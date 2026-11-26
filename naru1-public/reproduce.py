#!/usr/bin/env python3
"""Naru 1st-paper dry run: recompute the reviewer-strata tables from table.tsv alone.

Standard library only; reads nothing but the table next to this file (or argv[1]).
Columns: row, week, producer_tribe, producer, reviewer, reviewer_tribe,
first_decision, last_decision, n_decided, two_reviewers.
"""
import csv
import sys
from pathlib import Path


def load(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def approve(rows, key="last_decision"):
    return (sum(1 for r in rows if r[key] == "approve"), len(rows))


def summarize(rows):
    out = {}
    for t in sorted({r["producer_tribe"] for r in rows}):
        rs = [r for r in rows if r["producer_tribe"] == t]
        out[("A", t, "last")] = approve(rs)
        out[("A", t, "first")] = approve(rs, "first_decision")
    for k in sorted({(r["producer_tribe"], r["reviewer_tribe"]) for r in rows}):
        out[("B",) + k] = approve([r for r in rows if (r["producer_tribe"], r["reviewer_tribe"]) == k])
    for rv in sorted({r["reviewer"] for r in rows}):
        out[("D", rv)] = approve([r for r in rows if r["reviewer"] == rv])
    for k in sorted({(r["producer_tribe"], r["producer"], r["reviewer"]) for r in rows}):
        out[("H",) + k] = approve([r for r in rows if (r["producer_tribe"], r["producer"], r["reviewer"]) == k])
    out[("F", "two_reviewers")] = sum(1 for r in rows if r["two_reviewers"] == "yes")
    return out


def render(summary):
    lines = []
    for k, v in summary.items():
        head = "\t".join(k)
        if isinstance(v, tuple):
            lines.append(f"{head}\t{v[0]}/{v[1]}\t{100 * v[0] / v[1]:.1f}%")
        else:
            lines.append(f"{head}\t{v}")
    return "\n".join(lines)


if __name__ == "__main__":
    table = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("table.tsv")
    print(render(summarize(load(table))))
