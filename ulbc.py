from typing import List, Sequence, Tuple

# --------------------------------------------------------------------------- spec

SBOX = [0x8, 0x0, 0x1, 0x5, 0xC, 0x7, 0x4, 0x6,
        0x2, 0xA, 0x3, 0xD, 0xE, 0xF, 0xB, 0x9]

SBOX_INV = [0] * 16
for _i, _v in enumerate(SBOX):
    SBOX_INV[_v] = _i

# PostPerm:  Z[i] = Y[P_S[i]]
P_S = [0, 25, 18, 15, 4, 1, 30, 19, 8, 29, 22, 7, 12, 9, 26, 23,
       20, 5, 10, 31, 16, 13, 2, 27, 28, 17, 6, 3, 24, 21, 14, 11]

# Key schedule cell permutation:  K_i[j] = K_{i-1}[P_K[j]]
P_K = [8, 26, 4, 21, 9, 23, 15, 16, 20, 17, 11, 7, 24, 5, 25, 0,
       31, 18, 10, 27, 19, 2, 29, 13, 28, 6, 12, 22, 14, 30, 3, 1]

# uLBC-256 LFSR nibble substitution f1
F1 = [0x0, 0x8, 0x1, 0x9, 0x2, 0xA, 0x3, 0xB,
      0xC, 0x4, 0xD, 0x5, 0xE, 0x6, 0xF, 0x7]

ROUNDS = {128: 18, 256: 22}

NCELLS = 32
COLUMN_MAJOR = True


def columns(column_major: bool = COLUMN_MAJOR) -> List[Tuple[int, int, int, int]]:
    """The eight MixColumn cell-index groups under the chosen convention."""
    if column_major:
        return [(4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3) for j in range(8)]
    return [(j, j + 8, j + 16, j + 24) for j in range(8)]


COLS = columns(True)
COLS_ROWMAJOR = columns(False)

State = Tuple[int, ...]


# ----------------------------------------------------------------- round operations

def sub_nib(s: Sequence[int]) -> State:
    return tuple(SBOX[x] for x in s)


def sub_nib_inv(s: Sequence[int]) -> State:
    return tuple(SBOX_INV[x] for x in s)


def add_const(s: Sequence[int], rc: Sequence[int]) -> State:
    return tuple(x ^ c for x, c in zip(s, rc))


def post_perm(s: Sequence[int]) -> State:
    return tuple(s[P_S[i]] for i in range(NCELLS))


def post_perm_inv(s: Sequence[int]) -> State:
    out = [0] * NCELLS
    for i in range(NCELLS):
        out[P_S[i]] = s[i]
    return tuple(out)


def mix_column(s: Sequence[int], cols=COLS) -> State:
    out = list(s)
    for c in cols:
        tot = s[c[0]] ^ s[c[1]] ^ s[c[2]] ^ s[c[3]]
        for idx in c:
            out[idx] = tot ^ s[idx]
    return tuple(out)


mix_column_inv = mix_column  # involutory


def add_key(s: Sequence[int], k: Sequence[int]) -> State:
    return tuple(x ^ y for x, y in zip(s, k))


# ---------------------------------------------------------------------- key schedule

def perm_key(k: Sequence[int]) -> State:
    return tuple(k[P_K[j]] for j in range(NCELLS))


def key_schedule_128(mk: Sequence[int], rounds: int) -> List[State]:
    """Returns [K_{-1}, K_0, ..., K_{rounds-1}] (length rounds+1)."""
    keys = [tuple(mk)]
    cur = tuple(mk)
    for _ in range(rounds):
        cur = perm_key(cur)
        keys.append(cur)
    return keys


def key_schedule_256(mk0: Sequence[int], mk1: Sequence[int],
                     rounds: int) -> List[State]:
    out = []
    m1 = tuple(mk1)                             # f1^i(MK1), i = 0
    perm = tuple(mk0)                           # will hold F^i(.) incrementally
    for i in range(rounds + 1):
        combined = tuple(a ^ b for a, b in zip(mk0, m1))
        v = combined
        for _ in range(i):                      # F^i
            v = perm_key(v)
        out.append(v)
        m1 = tuple(F1[x] for x in m1)           # f1^{i+1}(MK1), nibble-wise
    del perm
    return out


# ------------------------------------------------------------------------- the cipher

class ULBC:
    def __init__(self, variant: int = 128, rounds: int = None,
                 rc: Sequence[Sequence[int]] = None, column_major: bool = True):
        self.variant = variant
        self.rounds = rounds if rounds is not None else ROUNDS[variant]
        self.cols = columns(column_major)
        self.rc = rc if rc is not None else [tuple([0] * NCELLS)] * (self.rounds + 1)

    # -- key material ------------------------------------------------------
    def expand(self, key) -> List[State]:
        if self.variant == 128:
            return key_schedule_128(key, self.rounds)
        mk0, mk1 = key
        return key_schedule_256(mk0, mk1, self.rounds)

    # -- single round ------------------------------------------------------
    def round_fwd(self, x: Sequence[int], k: Sequence[int], rc: Sequence[int],
                  last: bool = False) -> State:
        y = add_const(sub_nib(x), rc)
        if last:
            return add_key(y, k)
        z = post_perm(y)
        w = mix_column(z, self.cols)
        return add_key(w, k)

    def round_bwd(self, x: Sequence[int], k: Sequence[int], rc: Sequence[int],
                  last: bool = False) -> State:
        v = add_key(x, k)
        if last:
            return sub_nib_inv(add_const(v, rc))
        z = mix_column_inv(v, self.cols)
        y = post_perm_inv(z)
        return sub_nib_inv(add_const(y, rc))

    # -- full / reduced encryption ----------------------------------------
    def encrypt(self, p: Sequence[int], key, rounds: int = None,
                whiten: bool = True) -> State:
        r = self.rounds if rounds is None else rounds
        ks = self.expand(key)
        x = add_key(p, ks[0]) if whiten else tuple(p)
        for i in range(r):
            x = self.round_fwd(x, ks[i + 1], self.rc[i], last=(i == r - 1))
        return x

    def decrypt(self, c: Sequence[int], key, rounds: int = None,
                whiten: bool = True) -> State:
        r = self.rounds if rounds is None else rounds
        ks = self.expand(key)
        x = tuple(c)
        for i in range(r - 1, -1, -1):
            x = self.round_bwd(x, ks[i + 1], self.rc[i], last=(i == r - 1))
        return add_key(x, ks[0]) if whiten else x

    # -- trail access: the sequence X_0 .. X_r ----------------------------
    def trail(self, x0: Sequence[int], key, rounds: int) -> List[State]:
        ks = self.expand(key)
        xs = [tuple(x0)]
        x = tuple(x0)
        for i in range(rounds):
            x = self.round_fwd(x, ks[i + 1], self.rc[i], last=(i == rounds - 1))
            xs.append(x)
        return xs

    def trail_bwd(self, xr: Sequence[int], key, rounds: int,
                  start_round: int = 0) -> List[State]:
        """Invert `rounds` rounds ending at round index start_round+rounds-1.
        Returns states from X_{start_round} up to X_{start_round+rounds}."""
        ks = self.expand(key)
        xs = [tuple(xr)]
        x = tuple(xr)
        for i in range(start_round + rounds - 1, start_round - 1, -1):
            x = self.round_bwd(x, ks[i + 1], self.rc[i], last=False)
            xs.append(x)
        return list(reversed(xs))


# --------------------------------------------------------------------------- helpers

def rand_state(rng) -> State:
    return tuple(rng.randrange(16) for _ in range(NCELLS))


def active(delta: Sequence[int]) -> List[int]:
    return [i for i, d in enumerate(delta) if d]


def xor(a: Sequence[int], b: Sequence[int]) -> State:
    return tuple(x ^ y for x, y in zip(a, b))


def hexs(s: Sequence[int]) -> str:
    return ''.join('%x' % x for x in s)


def from_hex(h: str) -> State:
    return tuple(int(c, 16) for c in h)
