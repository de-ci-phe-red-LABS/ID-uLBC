import argparse
import json
from itertools import combinations

from ulbc import NCELLS, COLS, P_S
from exp1_distinguisher import (blank, cells_with, sym_round_fwd, sym_round_bwd,
                                sym_post_perm_inv, sym_mix_column, sym_sub_nib,
                                forward_trail, backward_trail, ZERO, NONZ, UNK)


# ------------------------------------------------------------------ pattern space
def single_cell_patterns():
    return [(c,) for c in range(NCELLS)]


def equal_pair_patterns():
    """Two cells in a common column carrying the same difference: the identical-
    difference trick the paper exploits, and the shape of its own output pattern."""
    out = []
    for col in COLS:
        for i, j in combinations(range(4), 2):
            out.append((col[i], col[j]))
    return out


def make_state(cells, symbol):
    st = blank()
    for c in cells:
        st[c] = (NONZ, frozenset({symbol}))
    return st


def fwd_states(cells, rounds):
    st = make_state(cells, 'a')
    out = [st]
    for _ in range(rounds):
        st = sym_round_fwd(st)
        out.append(st)
    return out


def bwd_states(cells, rounds):
    st = make_state(cells, 'b')
    out = [st]
    for _ in range(rounds):
        st = sym_round_bwd(st)
        out.append(st)
    return out          # out[k] = state k rounds back from the output


# ------------------------------------------------------------------ 2A  ID census
def find_ids(rounds, in_patterns, out_patterns, verbose=True):
    hits = []
    fwd_cache = {}
    bwd_cache = {}
    for ip in in_patterns:
        fwd_cache[ip] = fwd_states(ip, rounds - 1)
    for op in out_patterns:
        bwd_cache[op] = bwd_states(op, rounds - 1)

    for ip in in_patterns:
        f = fwd_cache[ip]
        for op in out_patterns:
            b = bwd_cache[op]
            for meet in range(1, rounds):
                fs = f[meet]
                bs = b[rounds - meet]
                contra = set(cells_with(fs, NONZ)) & set(cells_with(bs, ZERO))
                if contra:
                    hits.append({'in': list(ip), 'out': list(op), 'meet': meet,
                                 'cells': sorted(contra)})
                    break
    return hits


def census(max_rounds, verbose=True):
    singles = single_cell_patterns()
    pairs = equal_pair_patterns()
    table = {}
    print('\n2A  Census of deterministic impossible differentials')
    print('    input class : 1 active cell (%d patterns)' % len(singles))
    print('    output class: 1 active cell (%d) or equal co-columnar pair (%d)'
          % (len(singles), len(pairs)))
    print()
    print('    rounds   in=1,out=1        in=1,out=equal pair')
    print('    ------   ---------------   -----------------------')
    for r in range(4, max_rounds + 1):
        h11 = find_ids(r, singles, singles)
        h1p = find_ids(r, singles, pairs)
        table[r] = {'single_single': len(h11), 'single_pair': len(h1p),
                    'total_single_single': len(singles) ** 2,
                    'total_single_pair': len(singles) * len(pairs),
                    'examples_pair': h1p[:3], 'examples_single': h11[:3]}
        print('      %2d     %6d / %-6d   %6d / %-6d'
              % (r, len(h11), len(singles) ** 2, len(h1p), len(singles) * len(pairs)))
    return table


# ------------------------------------------------------------------ 2B  extensions
def extend_backward_one_round(cells):
    """Cells active in X_{-1} when X_0 has exactly `cells` active (deterministic
    part only).  Used to size the plaintext-side extension."""
    st = make_state(cells, 'x')
    back = sym_round_bwd(st)
    return (sorted(cells_with(back, NONZ)), sorted(cells_with(back, UNK)))


def extend_forward_one_round(cells):
    st = make_state(cells, 'y')
    fw = sym_round_fwd(st)
    return (sorted(cells_with(fw, NONZ)), sorted(cells_with(fw, UNK)))


def cost_extensions(dist, verbose=True):
    """For a distinguisher (in-pattern, out-pattern), report the boundary sizes of
    a 1-round backward and 1-round forward extension."""
    ip, op = tuple(dist['in']), tuple(dist['out'])
    bn, bu = extend_backward_one_round(ip)
    fn, fu = extend_forward_one_round(op)
    nb = len(bn) + len(bu)
    nf = len(fn) + len(fu)
    return {'in': list(ip), 'out': list(op), 'meet': dist['meet'],
            'plain_active': nb, 'plain_cells': sorted(bn + bu),
            'cipher_active': nf, 'cipher_cells': sorted(fn + fu),
            'delta_in_bits': 4 * nb, 'delta_out_bits': 4 * nf,
            # filtering: the backward extension must collapse nb active cells onto
            # the |ip| active cells of the distinguisher input, and likewise forward
            'cb_bits': 4 * (nb - len(ip)),
            'cf_bits': 4 * (nf - len(op))}


def rank_candidates(rounds, top=12):
    print('\n2B  Extension cost of the sparsest distinguishers at %d rounds' % rounds)
    singles = single_cell_patterns()
    pairs = equal_pair_patterns()
    hits = find_ids(rounds, singles, singles) + find_ids(rounds, singles, pairs)
    if not hits:
        print('    no distinguisher at %d rounds in this pattern class' % rounds)
        return []
    costed = [cost_extensions(h) for h in hits]
    costed.sort(key=lambda c: (c['cb_bits'] + c['cf_bits'],
                               c['delta_in_bits'] + c['delta_out_bits']))
    print('    in -> out                  meet  |Din|  |Dout|  c_b  c_f  c_b+c_f')
    print('    ------------------------   ----  -----  ------  ---  ---  -------')
    for c in costed[:top]:
        print('    %-24s   X_%-2d  %4d   %5d   %3d  %3d   %4d'
              % ('%s -> %s' % (c['in'], c['out']), c['meet'],
                 c['delta_in_bits'], c['delta_out_bits'],
                 c['cb_bits'], c['cf_bits'], c['cb_bits'] + c['cf_bits']))
    return costed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-rounds', type=int, default=10)
    ap.add_argument('--design-rounds', type=int, default=4,
                    help='round count of the distinguisher to cost for the '
                         'executable attack')
    ap.add_argument('--json', default='results_exp2.json')
    a = ap.parse_args()

    print('=' * 78)
    print('Stage 2: impossible-differential census and extension costing')
    print('=' * 78)

    table = census(a.max_rounds)

    longest = max((r for r, v in table.items()
                   if v['single_single'] or v['single_pair']), default=None)
    print('\n    longest distinguisher found in this pattern class: %s rounds'
          % longest)
    if longest is not None and longest >= 8:
        ex = table[longest]['examples_pair'] or table[longest]['examples_single']
        print('    example: input cell %s, output %s, contradiction at X_%d in cells %s'
              % (ex[0]['in'], ex[0]['out'], ex[0]['meet'], ex[0]['cells']))

    costed = {}
    for r in (a.design_rounds, a.design_rounds + 1):
        costed[r] = rank_candidates(r)

    with open(a.json, 'w') as fh:
        json.dump({'census': table, 'longest': longest,
                   'costed': {str(k): v for k, v in costed.items()}},
                  fh, indent=2, default=list)
    print('\nresults written to %s' % a.json)


if __name__ == '__main__':
    main()
