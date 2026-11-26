#!/usr/bin/env python3
"""Naru paper 1, v1 analysis: everything the manuscript computes from the frozen table alone.

Standard library only. Reads table.tsv next to this file (or argv[1]) and refuses any other bytes
than the frozen sha256. Prints the numbers the manuscript cites in sections 4.1-4.3.
"""
import csv
import hashlib
import sys
from collections import defaultdict
from math import comb
from pathlib import Path

FROZEN_SHA = "533bbeefe53af7581d8f7b8136c22b10efbc909a9a3d25082842269271c7c637"
path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().with_name("table.tsv")
data = path.read_bytes()
sha = hashlib.sha256(data).hexdigest()
if sha != FROZEN_SHA:
    sys.exit(f"refuse: table sha256 {sha} is not the frozen {FROZEN_SHA}")
rows = list(csv.DictReader(data.decode("utf-8").splitlines(), delimiter="\t"))


def ok(rs, key="last_decision"):
    return sum(r[key] == "approve" for r in rs), len(rs)


def pct(a, n):
    return f"{a}/{n} = {100 * a / n:.1f}%"


def fisher(a, b, c, d):
    r1, c1, n = a + b, a + c, a + b + c + d
    p = lambda x: comb(c1, x) * comb(n - c1, r1 - x) / comb(n, r1)  # noqa: E731
    obs = p(a)
    return sum(p(x) for x in range(max(0, r1 + c1 - n), min(r1, c1) + 1) if p(x) <= obs * (1 + 1e-9))


def wilson(a, n, z=1.959963984540054):
    ph = a / n
    den = 1 + z * z / n
    mid = (ph + z * z / (2 * n)) / den
    half = z * ((ph * (1 - ph) / n + z * z / (4 * n * n)) ** 0.5) / den
    return max(0.0, mid - half), min(1.0, mid + half)


print(f"table sha256 {sha} · rows {len(rows)}")
dec = defaultdict(int)
for r in rows:
    dec[r["last_decision"]] += 1
print("last decisions:", dict(sorted(dec.items())), "· rows with >1 decision:", sum(r["n_decided"] != "1" for r in rows))
print()
print("## 4.1 by producer tribe (pooled over reviewers) — the dashboard view")
for t in sorted({r["producer_tribe"] for r in rows}):
    rs = [r for r in rows if r["producer_tribe"] == t]
    print(f"  {t:11} last {pct(*ok(rs))}   first {pct(*ok(rs, 'first_decision'))}")
print()
print("## 4.2 by reviewer, and producer tribe within reviewer (Wilson 95% intervals)")
for rv in ("R1", "R2"):
    rs = [r for r in rows if r["reviewer"] == rv]
    a, n = ok(rs)
    lo, hi = wilson(a, n)
    print(f"  {rv} ({rs[0]['reviewer_tribe']}) all {pct(a, n)}  [{100 * lo:.1f}, {100 * hi:.1f}]")
    for t in sorted({r["producer_tribe"] for r in rs}):
        a, n = ok([r for r in rs if r["producer_tribe"] == t])
        lo, hi = wilson(a, n)
        print(f"     {t:11} {pct(a, n):18} [{100 * lo:.1f}, {100 * hi:.1f}]")
print()
print("## 4.2 reviewer contrast inside each tribe both reviewers saw (Fisher exact, two-sided)")
for t in ("codex_cli", "agy_cli"):
    a1, n1 = ok([r for r in rows if r["producer_tribe"] == t and r["reviewer"] == "R1"])
    a2, n2 = ok([r for r in rows if r["producer_tribe"] == t and r["reviewer"] == "R2"])
    print(f"  {t:11} R1 {a1}/{n1} vs R2 {a2}/{n2}  p = {fisher(a1, n1 - a1, a2, n2 - a2):.2e}")
print("  tribe contrast inside R1 (claude vs codex vs agy):",
      " · ".join(pct(*ok([r for r in rows if r['reviewer'] == 'R1' and r['producer_tribe'] == t]))
                 for t in ("claude_cli", "codex_cli", "agy_cli")))
a1, n1 = ok([r for r in rows if r["reviewer"] == "R1" and r["producer_tribe"] == "claude_cli"])
a2, n2 = ok([r for r in rows if r["reviewer"] == "R1" and r["producer_tribe"] == "codex_cli"])
print(f"  inside R1, claude {a1}/{n1} vs codex {a2}/{n2}  p = {fisher(a1, n1 - a1, a2, n2 - a2):.3f}")
a1, n1 = ok([r for r in rows if r["reviewer"] == "R2" and r["producer_tribe"] == "codex_cli"])
a2, n2 = ok([r for r in rows if r["reviewer"] == "R2" and r["producer_tribe"] == "agy_cli"])
print(f"  inside R2, codex {a1}/{n1} vs agy {a2}/{n2}  p = {fisher(a1, n1 - a1, a2, n2 - a2):.3f}")
print()
print("## 4.3 the empty cell and the pairs")
cells = defaultdict(int)
for r in rows:
    cells[(r["producer_tribe"], r["reviewer"])] += 1
print("  cells (producer tribe x reviewer):", dict(sorted(cells.items())))
print("  rows decided by both reviewers:", sum(r["two_reviewers"] == "yes" for r in rows))
print()
print("## 4.3 by ISO week of the deciding lease (reviewer x week)")
for wk in sorted({r["week"] for r in rows}):
    line = []
    for rv in ("R1", "R2"):
        rs = [r for r in rows if r["week"] == wk and r["reviewer"] == rv]
        line.append(f"{rv} {pct(*ok(rs)) if rs else '-':18}")
    print(f"  {wk}  " + "  ".join(line))
both = [wk for wk in sorted({r['week'] for r in rows})
        if {r["reviewer"] for r in rows if r["week"] == wk} == {"R1", "R2"}]
for rv in ("R1", "R2"):
    print(f"  weeks where both decided ({', '.join(both)}): {rv} {pct(*ok([r for r in rows if r['week'] in both and r['reviewer'] == rv]))}")
