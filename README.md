# ID-uLBC: Impossible Differential Cryptanalysis of uLBC

Companion code for the paper

> **Impossible Differential Attack on uLBC**
> Dilip Sau, Dhiman Saha, Dipanwita Roy Chowdhury
> INDOCRYPT 2026 (27th International Conference on Cryptology in India), Bengaluru, India, December 2026. Springer LNCS.

This repository complements the paper. It contains the constraint-programming model used to search for single-key impossible differential (ID) distinguishers of uLBC, the independent verification of the 8-round distinguisher reported in the paper, and the exhaustive census over the boundary class that reproduces the distinguisher outside the CP solver (Section 3 of the paper). The paper itself contains the analysis, the key-recovery attacks and all complexity estimates; nothing here is meant to be read without it.

## Contents

| File | Purpose | Paper section |
|---|---|---|
| `ulbc.py` | Reference implementation of uLBC (S-box, PostPerm, MixColumn, key schedules of uLBC-128 and uLBC-256). Used by the verification scripts. | 2.1 |
| `sk_distinguisher.mzn` | MiniZinc CP model of the uLBC round function in the cell-wise ZeroPlus style (activity pattern 0–3 plus difference value per nibble) for single-key ID search over `RD` rounds. | 3.1 |
| `run_all_sk.py` | Driver that instantiates the CP model for a range of round counts and records SAT/UNSAT and the trails found. | 3.1, 3.2 |
| `exp1_distinguisher.py` | Verification of the 8-round distinguisher (input nibble 31, equal output differences in nibbles 14 and 15): symbolic propagation, random-key empirical test, sweep over all 225 nonzero (a, b) pairs, and the 9-round check. | 3.2, 3.3 |
| `exp2_search.py` | Exhaustive census over the boundary class with one active input nibble and one or two (same-column, equal-difference) active output nibbles; reports the number of patterns giving an ID at each round count. | 3.3 |
| `results_exp1.json` | Output of `exp1_distinguisher.py` as used in the paper. | 3.3 |
| `results_exp2.json` | Output of `exp2_search.py` as used in the paper. | 3.3 |

## Requirements

* Python 3.9 or later. The verification scripts (`ulbc.py`, `exp1_distinguisher.py`, `exp2_search.py`) use only the standard library.
* For the CP search: [MiniZinc](https://www.minizinc.org/) 2.7 or later with the OR-Tools CP-SAT solver, and the `minizinc` Python package (`pip install minizinc`).

## Reproducing the results

**Verify the 8-round distinguisher** (symbolic + empirical, all 225 difference pairs, and the 9-round check):

```bash
python3 exp1_distinguisher.py --trials 20000 --sweep-trials 200 --json results_exp1.json
```

Expected outcome: forward trail forces nibbles {1, 3, 19, 21} of X4 to be active, backward trail forces nibbles {0, 1, 10, 11} of X4 to be inactive, contradiction in nibble 1; no failure over the 225 (a, b) pairs; 288 of the 1536 two-nibble output patterns yield an 8-round ID and none yields a 9-round ID.

**Run the exhaustive census** over the boundary class:

```bash
python3 exp2_search.py --max-rounds 10 --json results_exp2.json
```

Expected outcome: 288 of the 1536 single-input / equal-pair-output patterns give an 8-round ID, 0 of the 1024 single-input / single-output patterns do, and no pattern in the class gives a 9-round ID.

**Search with the CP model.** `run_all_sk.py` provides `test_sk(start_rd, end_rd)`, which runs `sk_distinguisher.mzn` for each round count in the range and writes the status and any trail found to `results_sk_all.txt`:

```bash
python3 -c "from run_all_sk import test_sk; test_sk(1, 10)"
```

The model is SAT for 8 rounds and UNSAT for 9 rounds under the input/output conditions described in the paper.

## Citation

If you use this code, please cite the paper:

```bibtex
@inproceedings{SauSahaRoyChowdhury2026uLBC,
  author    = {Dilip Sau and Dhiman Saha and Dipanwita {Roy Chowdhury}},
  title     = {Impossible Differential Attack on {uLBC}},
  booktitle = {Progress in Cryptology -- {INDOCRYPT} 2026},
  series    = {Lecture Notes in Computer Science},
  publisher = {Springer},
  year      = {2026}
}
```

## Authors

* Dilip Sau, Indian Institute of Technology Kharagpur
* Dhiman Saha, de.ci.phe.red Lab, Indian Institute of Technology Bhilai
* Dipanwita Roy Chowdhury, Indian Institute of Technology Kharagpur

## Acknowledgement

The CP model follows the cell-wise modelling approach of the ZeroPlus framework of Hadipour, Gerhalter, Sadeghi and Eichlseder (IACR ToSC 2024/1).
