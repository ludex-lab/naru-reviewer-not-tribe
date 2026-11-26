# The Reviewer, Not the Tribe — public bundle

Data and scripts for the paper *The Reviewer, Not the Tribe: Assignment Confounding in Cross-Model Audits Inside an AI Agent Society* (J. Jeong, November 2026).

The folder `naru1-public/` holds nine files:

| File | Role |
|---|---|
| `table.tsv` | The frozen, anonymized table: 231 decided landing audits |
| `reproduce.py` | Tables 1–2 and the within-tribe lines, from `table.tsv` alone |
| `analysis.py` | Wilson intervals, Fisher tests and the weekly split |
| `analysis_sensitivity.py` | The post hoc sensitivity analysis (not pre-registered), with a fixed seed |
| `reproduce.txt` · `analysis.txt` · `analysis_sensitivity.out` | Expected outputs |
| `SCHEMA.md` | Columns, removals, the checks run before freezing, what the schema cannot protect |
| `MANIFEST.sha256` | sha256 of the eight files above |

## Rerun

Python 3 standard library only. No network, no model calls.

```bash
cd naru1-public
shasum -a 256 -c MANIFEST.sha256
python3 -I -B reproduce.py | cmp - reproduce.txt && echo same
python3 -I -B analysis.py | cmp - analysis.txt && echo same
python3 -I -B analysis_sensitivity.py | cmp - analysis_sensitivity.out && echo same
```

Each script reads only the `table.tsv` beside it and refuses any other bytes. The sensitivity script takes about eight seconds.

## What is not here

The table was exported from a live request store by a script that reads unreleased files. The receipts behind the paper's sentences marked with a dagger (†) are not released. `SCHEMA.md` says what was removed and what the pseudonyms cannot protect.

## License

CC BY 4.0 (see `LICENSE`).
