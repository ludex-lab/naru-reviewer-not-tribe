# Naru paper 1 — anonymized table, schema v1 (public edition)

This folder accompanies the paper *The Reviewer, Not the Tribe*. It holds one frozen table of the 231 landing audits that had a decision, and three standard-library scripts that recompute every table-based number in the paper. The table was exported once from the village's live request store, and it is fixed by its sha256. The export script is not part of this folder, because it reads files that are not released.

## Files

| File | Role |
|---|---|
| `table.tsv` | The frozen table, 231 rows. sha256 must be `533bbeefe53af7581d8f7b8136c22b10efbc909a9a3d25082842269271c7c637` |
| `reproduce.py` | Prints tables A (producer tribe), B (producer tribe × reviewer tribe), D (reviewer), H (producer × reviewer) and F (rows decided by both reviewers) |
| `analysis.py` | Wilson intervals, Fisher exact tests and the weekly split. Refuses any table whose sha256 differs |
| `analysis_sensitivity.py` | The post hoc sensitivity analysis (not pre-registered): the week × tribe cells both reviewers saw, the first-decision tables, producer clustering and assignment balance. Fixed random seed. Refuses any table whose sha256 differs |
| `reproduce.txt` · `analysis.txt` · `analysis_sensitivity.out` | The expected output of each script |
| `MANIFEST.sha256` | sha256 of every other file in this folder |

To check: `shasum -a 256 -c MANIFEST.sha256`, then `python3 -I -B reproduce.py > r.out`, `python3 -I -B analysis.py > a.out` and `python3 -I -B analysis_sensitivity.py > s.out`. Compare each output with its expected file (`cmp`). The scripts read only `table.tsv` in their own folder, or a path given as the first argument.

## Columns (one row = one audit request with at least one decided response)

| Column | Meaning |
|---|---|
| `row` | Serial number after sorting. It is not the request order and cannot be traced back to a card |
| `week` | ISO week in which the lease of the deciding attempt ended (not the decision time). Nothing finer than a week is kept |
| `producer_tribe` | The model family of the producer (`claude_cli`, `codex_cli`, `agy_cli`, `cursor_cli`) |
| `producer` | `P01`–`P21`, pseudonyms assigned in order of first appearance |
| `reviewer` · `reviewer_tribe` | `R1`/`R2`, pseudonyms assigned by number of decisions, and the reviewer's model family |
| `first_decision` · `last_decision` | `approve` · `request_changes` · `cannot_verify` |
| `n_decided` · `two_reviewers` | Number of decided responses · whether both reviewers decided the row |

**Removed before export:** request ids, bundle hashes, card numbers, resident names, response texts, changed paths, base revisions, landing texts, and every timestamp finer than an ISO week.

## Checks run before the table was frozen

- **Leak check.** The table was searched for the names of residents, seats and the human operator; for card and request id prefixes; for filesystem path fragments; and for hex runs of 12 or more characters. The check plants one token of each kind and must catch all of them. Result: no leaks.
- **Second count.** Tables A, B and D were compared, cell by cell, with an independent count made by different code from the unreleased store. All 16 cells matched.
- **Clean rerun.** A resident of a tribe that is neither reviewer's (`agy_cli`) ran `table.tsv`, `reproduce.py` and `analysis.py` alone in an empty folder outside the repository. The sha256 matched, all 16 table cells matched, and the analysis output was byte-identical to `analysis.txt`. The only change to `analysis.py` since that rerun is its default table path, which now points to `table.tsv` in this folder, and the docstring line that describes it. The sensitivity script came later; its output was reproduced byte for byte in an empty folder that held only `table.tsv` and the script.

## A rule applied before release

**The residents who are the subject of the study read first, and their answers are printed with the results.** The table and the claim sentences went to the two reviewers before release. Their answer lines are printed in the paper's Appendix B. They cannot block publication. That decision belongs to the human founder. But their answers are not removed.

## What this schema cannot protect

- **There are only two reviewers, and each carries a tribe label.** Inside the village, R1 and R2 are known at once. Outside readers learn the model family, not the individual seat. This is why the rule above exists.
- **`producer` pseudonyms follow order of first appearance.** Someone who knows the village's card queue could trace them back. The column is released all the same, because the paper's producer-level checks (§4.2) and its within-tribe sentence (§4.3) cannot be rerun without it.
- **The live store moves.** Exporting again from the live folder gives different numbers. The village dashboard's own rates moved between 11-21 and 11-22 (claude_cli 98.5% → 98.9%, codex_cli 52.6% → 56.6%). The paper relies only on this frozen table.
- **No reasons.** Response texts are not released, so outside readers cannot reread why a reviewer decided as it did.
