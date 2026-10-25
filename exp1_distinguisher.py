import argparse
import json
import random
import sys
from collections import Counter

from ulbc import (ULBC, NCELLS, COLS, P_S, sub_nib, post_perm, mix_column,
                  rand_state, active, xor)

IN_CELL = 31
OUT_CELLS = (14, 15)

ZERO, NONZ, UNK = '0', 'N', '?'

# Symbolic deterministic difference propagation.

def sym_sub_nib(state):
    out = []
    for tag, form in state:
        if tag == ZERO:
            out.append((ZERO, None))
        elif tag == NONZ:
            out.append((NONZ, None))       # nonzero preserved, value unknown
        else:
            out.append((UNK, None))
    return out


def sym_post_perm(state):
    return [state[P_S[i]] for i in range(NCELLS)]


def sym_post_perm_inv(state):
    out = [None] * NCELLS
    for i in range(NCELLS):
        out[P_S[i]] = state[i]
    return out


def _xor_forms(cells):
    """XOR three cells symbolically.  Returns a cell tuple."""
    form = frozenset()
    for tag, f in cells:
        if tag == UNK:
            return (UNK, None)
        if tag == ZERO:
            continue
        if f is None:
            return (UNK, None)             # nonzero but unknown value: cannot combine
        form = form ^ f                    # XOR of linear forms = symmetric difference
    if not form:
        return (ZERO, None)                # deterministic cancellation
    return (NONZ, form)


def sym_mix_column(state):
    out = list(state)
    for col in COLS:
        for r in col:
            others = [state[i] for i in col if i != r]
            exact = _xor_forms(others)
            if exact[0] != UNK:
                out[r] = exact
                continue
            tags = [t for t, _ in others]
            if all(t == ZERO for t in tags):
                out[r] = (ZERO, None)
            elif tags.count(NONZ) == 1 and tags.count(UNK) == 0:
                out[r] = (NONZ, None)
            else:
                out[r] = (UNK, None)
    return out


def sym_round_fwd(state):
    return sym_mix_column(sym_post_perm(sym_sub_nib(state)))


def sym_round_bwd(state):
    # inverse of (SB, AC, PP, MC, AK): AK and AC vanish on differences
    return sym_sub_nib(sym_post_perm_inv(sym_mix_column(state)))


def blank():
    return [(ZERO, None)] * NCELLS


def tags(state):
    return ''.join(t for t, _ in state)


def cells_with(state, tag):
    return [i for i, (t, _) in enumerate(state) if t == tag]


def forward_trail(in_cell=IN_CELL, rounds=4):
    st = blank()
    st[in_cell] = (NONZ, frozenset({'a'}))
    out = [st]
    for _ in range(rounds):
        st = sym_round_fwd(st)
        out.append(st)
    return out


def backward_trail(out_cells=OUT_CELLS, rounds=4, equal=True):
    st = blank()
    for k, c in enumerate(out_cells):
        st[c] = (NONZ, frozenset({'b'}) if equal else frozenset({'b%d' % k}))
    out = [st]
    for _ in range(rounds):
        st = sym_round_bwd(st)
        out.append(st)
    return list(reversed(out))


def run_symbolic(rd=8, meet=4, verbose=True):
    fwd = forward_trail(IN_CELL, meet)
    bwd = backward_trail(OUT_CELLS, rd - meet)

    if verbose:
        print('\n1A  Symbolic deterministic propagation (independent of ZeroPlus)')
        print('    forward from dX_0[%d] = a' % IN_CELL)
        for r, st in enumerate(fwd):
            print('      X_%d  certainly-nonzero: %2d   certainly-zero: %2d   unknown: %2d'
                  % (r, len(cells_with(st, NONZ)), len(cells_with(st, ZERO)),
                     len(cells_with(st, UNK))))
        print('    backward from dX_%d[%d] = dX_%d[%d] = b' %
              (rd, OUT_CELLS[0], rd, OUT_CELLS[1]))
        for k, st in enumerate(bwd):
            r = meet + k
            print('      X_%d  certainly-nonzero: %2d   certainly-zero: %2d   unknown: %2d'
                  % (r, len(cells_with(st, NONZ)), len(cells_with(st, ZERO)),
                     len(cells_with(st, UNK))))

    f_nonzero = set(cells_with(fwd[meet], NONZ))
    b_zero = set(cells_with(bwd[0], ZERO))
    contradiction = sorted(f_nonzero & b_zero)

    if verbose:
        print('    at X_%d:' % meet)
        print('      forward says certainly NONZERO : %s' % sorted(f_nonzero))
        print('      backward says certainly ZERO   : %s' % sorted(b_zero))
        print('      CONTRADICTION cells            : %s' % contradiction)
    return {'forward_nonzero': sorted(f_nonzero), 'backward_zero': sorted(b_zero),
            'contradiction': contradiction,
            'fwd_counts': [len(cells_with(s, NONZ)) for s in fwd],
            'bwd_counts': [len(cells_with(s, NONZ)) for s in bwd]}


# =====================================================================  1B
def full_rounds_fwd(cipher, ks, x, r0, nr):
    """Apply rounds r0 .. r0+nr-1 as full rounds (never the truncated last round)."""
    for i in range(r0, r0 + nr):
        x = cipher.round_fwd(x, ks[i + 1], cipher.rc[i], last=False)
    return x


def full_rounds_bwd(cipher, ks, x, r0, nr):
    """Invert rounds r0+nr-1 .. r0."""
    for i in range(r0 + nr - 1, r0 - 1, -1):
        x = cipher.round_bwd(x, ks[i + 1], cipher.rc[i], last=False)
    return x


def run_empirical(trials, meet=4, rd=8, seed=1, verbose=True):
    rng = random.Random(seed)
    cipher = ULBC(128)

    fwd_active_count = Counter()
    bwd_active_count = Counter()
    fwd_counts_per_round = [Counter() for _ in range(meet + 1)]

    for _ in range(trials):
        mk = rand_state(rng)
        ks = cipher.expand(mk)

        # ---- forward: only cell IN_CELL differs in X_0
        x = rand_state(rng)
        a = rng.randrange(1, 16)
        xp = list(x)
        xp[IN_CELL] ^= a
        xp = tuple(xp)
        u, v = x, xp
        for r in range(meet):
            u = cipher.round_fwd(u, ks[r + 1], cipher.rc[r], last=False)
            v = cipher.round_fwd(v, ks[r + 1], cipher.rc[r], last=False)
            d = xor(u, v)
            fwd_counts_per_round[r + 1][len(active(d))] += 1
        for c in active(xor(u, v)):
            fwd_active_count[c] += 1

        # ---- backward: cells OUT_CELLS carry the same difference b in X_rd
        y = rand_state(rng)
        b = rng.randrange(1, 16)
        yp = list(y)
        for c in OUT_CELLS:
            yp[c] ^= b
        yp = tuple(yp)
        u2, v2 = y, tuple(yp)
        u2 = full_rounds_bwd(cipher, ks, u2, meet, rd - meet)
        v2 = full_rounds_bwd(cipher, ks, v2, meet, rd - meet)
        for c in active(xor(u2, v2)):
            bwd_active_count[c] += 1

    always_active_fwd = sorted(c for c in range(NCELLS)
                               if fwd_active_count[c] == trials)
    never_active_bwd = sorted(c for c in range(NCELLS) if bwd_active_count[c] == 0)
    contradiction = sorted(set(always_active_fwd) & set(never_active_bwd))

    if verbose:
        print('\n1B  Empirical verification on the real cipher (%d random trials)' % trials)
        print('    active-cell counts observed forward from a 1-cell difference:')
        for r in range(1, meet + 1):
            dist = fwd_counts_per_round[r]
            lo, hi = min(dist), max(dist)
            print('      X_%d : %s' % (r, ('always %d' % lo) if lo == hi
                                        else 'range %d..%d' % (lo, hi)))
        print('    cells ALWAYS active at X_%d (forward)  : %s'
              % (meet, always_active_fwd))
        print('    cells NEVER  active at X_%d (backward) : %s'
              % (meet, never_active_bwd))
        print('    experimentally contradicting cells     : %s' % contradiction)
    return {'always_active_fwd': always_active_fwd,
            'never_active_bwd': never_active_bwd,
            'contradiction': contradiction,
            'trials': trials}


# =====================================================================  1C
def run_sweep(trials_per_pair, meet=4, rd=8, seed=7, verbose=True):
    rng = random.Random(seed)
    cipher = ULBC(128)
    worst = None
    failures = []

    for a in range(1, 16):
        for b in range(1, 16):
            fwd_always = set(range(NCELLS))
            bwd_ever = set()
            for _ in range(trials_per_pair):
                mk = rand_state(rng)
                ks = cipher.expand(mk)

                x = rand_state(rng)
                xp = list(x)
                xp[IN_CELL] ^= a
                u = full_rounds_fwd(cipher, ks, x, 0, meet)
                v = full_rounds_fwd(cipher, ks, tuple(xp), 0, meet)
                fwd_always &= set(active(xor(u, v)))

                y = rand_state(rng)
                yp = list(y)
                for c in OUT_CELLS:
                    yp[c] ^= b
                u2 = full_rounds_bwd(cipher, ks, y, meet, rd - meet)
                v2 = full_rounds_bwd(cipher, ks, tuple(yp), meet, rd - meet)
                bwd_ever |= set(active(xor(u2, v2)))

            contra = sorted(fwd_always - bwd_ever)
            if not contra:
                failures.append((a, b))
            if worst is None or len(contra) < worst[2]:
                worst = (a, b, len(contra), contra)

    if verbose:
        print('\n1C  Sweep over all 225 nonzero (a, b) pairs, %d trials each'
              % trials_per_pair)
        print('    pairs with NO contradiction (would break the distinguisher): %d'
              % len(failures))
        print('    weakest pair: a=%x b=%x with %d contradicting cell(s): %s'
              % (worst[0], worst[1], worst[2], worst[3]))
    return {'failures': failures, 'weakest': list(worst[:3]) + [worst[3]],
            'trials_per_pair': trials_per_pair}


# =====================================================================  1D
def run_nine_round(verbose=True):
    if verbose:
        print('\n1D  The 9-round question (within this deterministic model)')

    res9 = run_symbolic(rd=9, meet=4, verbose=False)
    if verbose:
        print('    same input/output, 9 rounds (meet at X_4): contradiction cells %s'
              % res9['contradiction'])

    # sweep meeting points for 9 rounds with the paper's boundary conditions
    best = []
    for meet in range(1, 9):
        fwd = forward_trail(IN_CELL, meet)
        bwd = backward_trail(OUT_CELLS, 9 - meet)
        contra = sorted(set(cells_with(fwd[meet], NONZ)) &
                        set(cells_with(bwd[0], ZERO)))
        best.append((meet, contra))
    if verbose:
        print('    9 rounds, contradiction by meeting point:')
        for meet, contra in best:
            print('      meet at X_%d : %s' % (meet, contra if contra else '(none)'))

    # broader: every single active input cell x every equal-difference output pair
    # in a common column, for 8 and 9 rounds
    found = {8: [], 9: []}
    col_pairs = [(c[i], c[j]) for c in COLS for i in range(4) for j in range(i + 1, 4)]
    for rd in (8, 9):
        for ic in range(NCELLS):
            for oc in col_pairs:
                for meet in range(1, rd):
                    fwd = forward_trail(ic, meet)
                    bwd = backward_trail(oc, rd - meet)
                    if set(cells_with(fwd[meet], NONZ)) & set(cells_with(bwd[0], ZERO)):
                        found[rd].append((ic, oc, meet))
                        break
    if verbose:
        print('    exhaustive search over 32 input cells x %d co-columnar equal output'
              % len(col_pairs))
        print('    pairs x all meeting points:')
        for rd in (8, 9):
            print('      %d rounds: %d of %d boundary patterns give a deterministic '
                  'contradiction' % (rd, len(found[rd]), NCELLS * len(col_pairs)))
    return {'nine_same_pattern': res9['contradiction'],
            'nine_by_meet': [(m, c) for m, c in best],
            'patterns_8': len(found[8]), 'patterns_9': len(found[9]),
            'total_patterns': NCELLS * len(col_pairs)}


# =====================================================================  main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--trials', type=int, default=20000)
    ap.add_argument('--sweep-trials', type=int, default=200)
    ap.add_argument('--json', default='results_exp1.json')
    a = ap.parse_args()

    print('=' * 78)
    print('Stage 1: validation of the 8-round impossible differential of uLBC')
    print('=' * 78)
    print('input   : dX_0[%d] = a != 0' % IN_CELL)
    print('output  : dX_8[%d] = dX_8[%d] = b != 0' % OUT_CELLS)

    sym = run_symbolic()
    emp = run_empirical(a.trials)
    swp = run_sweep(a.sweep_trials)
    nine = run_nine_round()

    print('\n' + '=' * 78)
    print('VERDICT')
    agree = set(sym['contradiction']) & set(emp['contradiction'])
    if sym['contradiction'] and emp['contradiction'] and agree:
        print('  The 8-round differential is IMPOSSIBLE.')
        print('  Symbolic and empirical analyses agree on contradiction cell(s): %s'
              % sorted(agree))
        print('  Holds for all 225 nonzero (a,b) pairs: %s'
              % ('yes' if not swp['failures'] else 'NO -- see failures'))
    else:
        print('  COULD NOT CONFIRM the distinguisher. Symbolic: %s  Empirical: %s'
              % (sym['contradiction'], emp['contradiction']))
    print('=' * 78)

    with open(a.json, 'w') as fh:
        json.dump({'symbolic': sym, 'empirical': emp, 'sweep': swp, 'nine': nine},
                  fh, indent=2, default=list)
    print('results written to %s' % a.json)


if __name__ == '__main__':
    main()
