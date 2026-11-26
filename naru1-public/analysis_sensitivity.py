#!/usr/bin/env python3
"""Post-hoc sensitivity analysis for Paper 1 ("The Reviewer, Not the Tribe").

Strictly reads ONLY the SHA-256-frozen public table
``table.tsv`` next to this file (231 rows).
Python standard library only. Zero LLM calls. Deterministic fixed RNG seed.
Uses anonymized IDs only (R1, R2, P01-P21); no resident names, paths, or cards.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import math
from pathlib import Path
import random
import sys

FROZEN_SHA = "533bbeefe53af7581d8f7b8136c22b10efbc909a9a3d25082842269271c7c637"
WILSON_Z = 1.959963984540054
BOOT_B = 20000
PERM_B = 100000
SEED = 20261125


def _default_table_path() -> Path:
    local = Path(__file__).resolve().with_name("table.tsv")
    if local.is_file():
        return local
    sys.exit("refuse: no table.tsv next to this file; pass the table path as the first argument")


def load_frozen_table(path: Path) -> tuple[str, list[dict[str, str]]]:
    data = path.read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    if sha != FROZEN_SHA:
        sys.exit(f"refuse: table sha256 {sha} is not the frozen {FROZEN_SHA}")
    rows = list(csv.DictReader(data.decode("utf-8").splitlines(), delimiter="\t"))
    return sha, rows


def wilson(a: int, n: int, z: float = WILSON_Z) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    ph = a / n
    den = 1.0 + z * z / n
    mid = (ph + z * z / (2.0 * n)) / den
    half = z * math.sqrt(ph * (1.0 - ph) / n + z * z / (4.0 * n * n)) / den
    return max(0.0, mid - half), min(1.0, mid + half)


def fisher_two_sided(a: int, b: int, c: int, d: int) -> float:
    r1, c1, n = a + b, a + c, a + b + c + d
    if n == 0:
        return 1.0
    denom = math.comb(n, r1)
    lo = max(0, r1 + c1 - n)
    hi = min(r1, c1)
    probs = {x: math.comb(c1, x) * math.comb(n - c1, r1 - x) / denom for x in range(lo, hi + 1)}
    obs = probs[a]
    return min(1.0, sum(p for p in probs.values() if p <= obs * (1.0 + 1e-12)))


def _convolve_dists(dists: list[dict[int, float]]) -> dict[int, float]:
    acc: dict[int, float] = {0: 1.0}
    for d in dists:
        nxt: dict[int, float] = defaultdict(float)
        for s0, p0 in acc.items():
            for x, px in d.items():
                nxt[s0 + x] += p0 * px
        acc = dict(nxt)
    return acc


def _noncentral_hypergeom_dist(n1: int, n2: int, m1: int, log_or: float) -> dict[int, float]:
    lo = max(0, m1 - n2)
    hi = min(n1, m1)
    log_w: dict[int, float] = {}
    for x in range(lo, hi + 1):
        lw = math.lgamma(n1 + 1) - math.lgamma(x + 1) - math.lgamma(n1 - x + 1)
        lw += math.lgamma(n2 + 1) - math.lgamma(m1 - x + 1) - math.lgamma(n2 - (m1 - x) + 1)
        lw += x * log_or
        log_w[x] = lw
    max_lw = max(log_w.values())
    w = {x: math.exp(lw - max_lw) for x, lw in log_w.items()}
    total = sum(w.values())
    return {x: val / total for x, val in w.items()}


def exact_cmh(strata_tables: list[tuple[int, int, int, int]]) -> dict[str, float]:
    """Exact conditional CMH test, MH common OR (RBG 95% CI), and exact conditional OR + 95% CI."""
    # Mantel-Haenszel OR + Robins-Breslow-Greenland (1986) 95% CI
    num_mh = 0.0
    den_mh = 0.0
    rd_num = 0.0
    rd_den = 0.0
    sum_pr = 0.0
    sum_ps_qr = 0.0
    sum_qs = 0.0
    s_obs = 0
    s_min = 0
    s_max = 0

    null_dists: list[dict[int, float]] = []
    for a, b, c, d in strata_tables:
        n1 = a + b
        n2 = c + d
        m1 = a + c
        ntot = n1 + n2
        if ntot == 0:
            continue
        r_k = a * d / ntot
        s_k = b * c / ntot
        p_k = (a + d) / ntot
        q_k = (b + c) / ntot
        num_mh += r_k
        den_mh += s_k
        sum_pr += p_k * r_k
        sum_ps_qr += p_k * s_k + q_k * r_k
        sum_qs += q_k * s_k
        w_rd = (n1 * n2) / ntot
        rd_num += w_rd * ((a / n1) - (c / n2))
        rd_den += w_rd
        s_obs += a
        s_min += max(0, m1 - n2)
        s_max += min(n1, m1)
        null_dists.append(_noncentral_hypergeom_dist(n1, n2, m1, 0.0))

    or_mh = num_mh / den_mh if den_mh > 0 else float("inf")
    if num_mh > 0 and den_mh > 0:
        var_log_or = (
            sum_pr / (2.0 * num_mh * num_mh)
            + sum_ps_qr / (2.0 * num_mh * den_mh)
            + sum_qs / (2.0 * den_mh * den_mh)
        )
        se_log_or = math.sqrt(var_log_or)
        rbg_lo = math.exp(math.log(or_mh) - WILSON_Z * se_log_or)
        rbg_hi = math.exp(math.log(or_mh) + WILSON_Z * se_log_or)
    else:
        rbg_lo, rbg_hi = float("nan"), float("nan")

    rd_mh = rd_num / rd_den if rd_den > 0 else 0.0

    # Exact conditional null distribution of S = sum_k a_k
    conv_null = _convolve_dists(null_dists)
    p_obs = conv_null[s_obs]
    p_two_sided = min(1.0, sum(p for p in conv_null.values() if p <= p_obs * (1.0 + 1e-12)))
    p_upper = sum(p for s, p in conv_null.items() if s >= s_obs)

    def _mean_s(log_or: float) -> float:
        m = 0.0
        for a, b, c, d in strata_tables:
            dist = _noncentral_hypergeom_dist(a + b, c + d, a + c, log_or)
            m += sum(x * p for x, p in dist.items())
        return m

    def _tail_ge(log_or: float) -> float:
        dists = [_noncentral_hypergeom_dist(a + b, c + d, a + c, log_or) for a, b, c, d in strata_tables]
        conv = _convolve_dists(dists)
        return sum(p for s, p in conv.items() if s >= s_obs)

    def _tail_le(log_or: float) -> float:
        dists = [_noncentral_hypergeom_dist(a + b, c + d, a + c, log_or) for a, b, c, d in strata_tables]
        conv = _convolve_dists(dists)
        return sum(p for s, p in conv.items() if s <= s_obs)

    def _bisect(fn, target: float, lo: float = -20.0, hi: float = 20.0) -> float:
        f_lo = fn(lo) - target
        for _ in range(90):
            mid = 0.5 * (lo + hi)
            f_mid = fn(mid) - target
            if (f_lo <= 0 and f_mid <= 0) or (f_lo > 0 and f_mid > 0):
                lo = mid
                f_lo = f_mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    if s_obs == s_max:
        cmle_or = float("inf")
        exact_hi = float("inf")
    elif s_obs == s_min:
        cmle_or = 0.0
        exact_lo = 0.0
        exact_hi = math.exp(_bisect(_tail_le, 0.025))
    else:
        cmle_or = math.exp(_bisect(_mean_s, float(s_obs)))
        exact_hi = math.exp(_bisect(_tail_le, 0.025))
    if s_obs > s_min:
        exact_lo = math.exp(_bisect(_tail_ge, 0.025))

    return {
        "s_obs": float(s_obs),
        "s_min": float(s_min),
        "s_max": float(s_max),
        "or_mh": or_mh,
        "rbg_lo": rbg_lo,
        "rbg_hi": rbg_hi,
        "cmle_or": cmle_or,
        "exact_lo": exact_lo,
        "exact_hi": exact_hi,
        "p_two_sided": p_two_sided,
        "p_upper": p_upper,
        "rd_mh": rd_mh,
    }


def _pct(a: int, n: int) -> str:
    return f"{a}/{n} = {100.0 * a / n:.1f}%" if n > 0 else "0/0 = n/a"


def _fmt_or(val: float) -> str:
    if math.isinf(val):
        return "inf"
    if math.isnan(val):
        return "nan"
    return f"{val:.2f}"


def _chi2_homogeneity(counts: list[tuple[int, int]]) -> float:
    tot_r1 = sum(x for x, _ in counts)
    tot_r2 = sum(y for _, y in counts)
    n = tot_r1 + tot_r2
    if n == 0 or tot_r1 == 0 or tot_r2 == 0:
        return 0.0
    p1 = tot_r1 / n
    p2 = tot_r2 / n
    chi2 = 0.0
    for r1, r2 in counts:
        m = r1 + r2
        if m == 0:
            continue
        e1 = m * p1
        e2 = m * p2
        chi2 += (r1 - e1) ** 2 / e1 + (r2 - e2) ** 2 / e2
    return chi2


def perm_homogeneity_test(counts: list[tuple[int, int]], *, n_perm: int, rng: random.Random) -> tuple[float, float]:
    """Monte Carlo exact conditional permutation test for K x 2 assignment homogeneity."""
    obs = _chi2_homogeneity(counts)
    row_sizes = [r1 + r2 for r1, r2 in counts]
    tot_r1 = sum(r1 for r1, _ in counts)
    n = sum(row_sizes)
    if len(counts) <= 1 or tot_r1 == 0 or tot_r1 == n:
        return obs, 1.0
    p1 = tot_r1 / n
    p2 = 1.0 - p1
    labels = [1] * tot_r1 + [0] * (n - tot_r1)
    ge = 0
    for _ in range(n_perm):
        rng.shuffle(labels)
        pos = 0
        stat = 0.0
        for m in row_sizes:
            r1 = sum(labels[pos : pos + m])
            pos += m
            r2 = m - r1
            e1 = m * p1
            e2 = m * p2
            stat += (r1 - e1) ** 2 / e1 + (r2 - e2) ** 2 / e2
        if stat >= obs - 1e-12:
            ge += 1
    return obs, (ge + 1) / (n_perm + 1)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    path = Path(args[0]) if args else _default_table_path()
    sha, rows = load_frozen_table(path)

    print("사후 민감도 분석, 사전등록 아님 (Post-hoc sensitivity analysis; not pre-registered)")
    print(f"source_table_sha256: {sha} · rows: {len(rows)} · rng_seed: {SEED}")
    print("anonymization_check: IDs limited to R1/R2 and P01-P21; zero resident names, paths, or card IDs.")
    print()

    # ------------------------------------------------------------------
    # (a) Overlapping strata (week x producer_tribe) — Review §3.1
    # ------------------------------------------------------------------
    print("## (a) 겹치는 칸만의 층화 비교 — (week x producer_tribe) [Review §3.1]")

    def _run_section_a(subset: list[dict[str, str]], label: str) -> None:
        by_cell: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for r in subset:
            by_cell[(r["week"], r["producer_tribe"])].append(r)
        overlap_keys = sorted(
            k for k, rs in by_cell.items()
            if {r["reviewer"] for r in rs} == {"R1", "R2"}
        )
        print(f"### (a) {label} (overlapping strata = {len(overlap_keys)}, rows in overlap = {sum(len(by_cell[k]) for k in overlap_keys)})")
        tables: list[tuple[int, int, int, int]] = []
        tot_r1_ok = tot_r1_n = tot_r2_ok = tot_r2_n = 0
        for wk, trb in overlap_keys:
            rs = by_cell[(wk, trb)]
            r1_rs = [r for r in rs if r["reviewer"] == "R1"]
            r2_rs = [r for r in rs if r["reviewer"] == "R2"]
            a = sum(r["last_decision"] == "approve" for r in r1_rs)
            b = len(r1_rs) - a
            c = sum(r["last_decision"] == "approve" for r in r2_rs)
            d = len(r2_rs) - c
            tables.append((a, b, c, d))
            tot_r1_ok += a
            tot_r1_n += len(r1_rs)
            tot_r2_ok += c
            tot_r2_n += len(r2_rs)
            p_cell = fisher_two_sided(a, b, c, d)
            diff_pp = 100.0 * (a / len(r1_rs) - c / len(r2_rs))
            print(f"  {wk}  {trb:10}  R1 {_pct(a, len(r1_rs)):18}  R2 {_pct(c, len(r2_rs)):18}  diff={diff_pp:+5.1f} pp  Fisher p={p_cell:.4g}")
        res = exact_cmh(tables)
        print(f"  Pooled over overlap strata: R1 {_pct(tot_r1_ok, tot_r1_n)} vs R2 {_pct(tot_r2_ok, tot_r2_n)}")
        print(
            f"  Mantel-Haenszel OR = {_fmt_or(res['or_mh'])} "
            f"(RBG 95% CI [{_fmt_or(res['rbg_lo'])}, {_fmt_or(res['rbg_hi'])}]) · "
            f"MH weighted RD = {100.0 * res['rd_mh']:+.1f} pp"
        )
        print(
            f"  Exact conditional CMH: S_obs={int(res['s_obs'])} (support [{int(res['s_min'])}, {int(res['s_max'])}]) · "
            f"CMLE OR = {_fmt_or(res['cmle_or'])} (exact 95% CI [{_fmt_or(res['exact_lo'])}, {_fmt_or(res['exact_hi'])}]) · "
            f"exact two-sided p = {res['p_two_sided']:.4e} (one-sided upper p = {res['p_upper']:.4e})"
        )

    _run_section_a(rows, "전체 결정 행 (All decisions)")
    no_cv = [r for r in rows if r["last_decision"] != "cannot_verify"]
    print(f"  [note] rows removed as cannot_verify across entire table: {len(rows) - len(no_cv)} (remaining {len(no_cv)})")
    _run_section_a(no_cv, "cannot_verify 제외 (Excluding cannot_verify)")
    print()

    # ------------------------------------------------------------------
    # (b) Decision-timing sensitivity — Review §3.4
    # ------------------------------------------------------------------
    print("## (b) 결정 시점 민감도 — 표 2 재집계 (Wilson 95% CI) [Review §3.4]")

    def _print_table2(spec_rows: list[tuple[str, str, str]], title: str) -> None:
        print(f"### {title} (n = {len(spec_rows)})")
        for rv in ("R1", "R2"):
            rs = [x for x in spec_rows if x[0] == rv]
            a = sum(x[2] == "approve" for x in rs)
            lo, hi = wilson(a, len(rs))
            print(f"  {rv} all        {_pct(a, len(rs)):18}  Wilson 95% CI [{100 * lo:.1f}, {100 * hi:.1f}]")
            for trb in sorted({x[1] for x in rs}):
                sub = [x for x in rs if x[1] == trb]
                sa = sum(x[2] == "approve" for x in sub)
                slo, shi = wilson(sa, len(sub))
                print(f"     {trb:10} {_pct(sa, len(sub)):18}  Wilson 95% CI [{100 * slo:.1f}, {100 * shi:.1f}]")
        for trb in ("codex_cli", "agy_cli"):
            r1s = [x for x in spec_rows if x[0] == "R1" and x[1] == trb]
            r2s = [x for x in spec_rows if x[0] == "R2" and x[1] == trb]
            a1 = sum(x[2] == "approve" for x in r1s)
            a2 = sum(x[2] == "approve" for x in r2s)
            p_f = fisher_two_sided(a1, len(r1s) - a1, a2, len(r2s) - a2)
            print(f"  contrast {trb:10}: R1 {a1}/{len(r1s)} vs R2 {a2}/{len(r2s)}  Fisher two-sided p = {p_f:.4e}")

    single_dec_rows = [r for r in rows if r["n_decided"] == "1"]
    assert all(r["first_decision"] == r["last_decision"] for r in single_dec_rows)
    _print_table2(
        [(r["reviewer"], r["producer_tribe"], r["last_decision"]) for r in single_dec_rows],
        "(b1) n_decided == 1 행만 (Single-decision requests only)",
    )

    two_rv_rows = [r for r in rows if r["two_reviewers"] == "yes"]
    print(
        f"  [note] rows with two_reviewers == 'yes': {len(two_rv_rows)} -> "
        + ", ".join(
            f"row {r['row']} ({r['week']}, {r['producer_tribe']}, {r['producer']}, "
            f"last_reviewer={r['reviewer']}, first={r['first_decision']}, last={r['last_decision']})"
            for r in two_rv_rows
        )
    )
    first_spec: list[tuple[str, str, str]] = []
    for r in rows:
        rv = r["reviewer"]
        if r["two_reviewers"] == "yes":
            rv = "R1" if rv == "R2" else "R2"
        first_spec.append((rv, r["producer_tribe"], r["first_decision"]))
    _print_table2(
        first_spec,
        "(b2) 첫 결정 기준 (First decision; two_reviewers='yes' row attributed to first reviewer)",
    )
    print()

    # ------------------------------------------------------------------
    # (c) Producer clustering — Review §3.3
    # ------------------------------------------------------------------
    print("## (c) 생산자 군집 민감도 — 생산자 내 짝비교 및 Cluster Bootstrap [Review §3.3]")
    by_prod: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_prod[r["producer"]].append(r)

    both_prods = sorted(
        p for p, rs in by_prod.items()
        if {r["reviewer"] for r in rs} == {"R1", "R2"}
    )
    print(f"### (c1) 두 검토자를 모두 만난 생산자 내 R1 - R2 승인율 차이 (n_producers = {len(both_prods)})")
    diffs: list[float] = []
    for p in both_prods:
        rs = by_prod[p]
        trb = rs[0]["producer_tribe"]
        r1s = [r for r in rs if r["reviewer"] == "R1"]
        r2s = [r for r in rs if r["reviewer"] == "R2"]
        a1 = sum(r["last_decision"] == "approve" for r in r1s)
        a2 = sum(r["last_decision"] == "approve" for r in r2s)
        d_val = (a1 / len(r1s)) - (a2 / len(r2s))
        diffs.append(d_val)
        print(f"  {p} ({trb:10}): R1 {_pct(a1, len(r1s)):17}  R2 {_pct(a2, len(r2s)):17}  R1-R2 = {100.0 * d_val:+6.1f} pp")

    s_diffs = sorted(diffs)
    n_p = len(s_diffs)
    med = s_diffs[n_p // 2] if n_p % 2 == 1 else 0.5 * (s_diffs[n_p // 2 - 1] + s_diffs[n_p // 2])
    mean_d = sum(s_diffs) / n_p
    pos = sum(x > 0 for x in s_diffs)
    zero = sum(x == 0 for x in s_diffs)
    neg = sum(x < 0 for x in s_diffs)
    m_nz = pos + neg
    k_small = min(pos, neg)
    p_sign = min(1.0, 2.0 * sum(math.comb(m_nz, k) for k in range(k_small + 1)) / (2 ** m_nz)) if m_nz > 0 else 1.0
    print(
        f"  Paired diff summary (n={n_p}): pos={pos}, zero={zero}, neg={neg} · "
        f"min={100*s_diffs[0]:+.1f} pp, median={100*med:+.1f} pp, mean={100*mean_d:+.1f} pp, max={100*s_diffs[-1]:+.1f} pp"
    )
    print(f"  Exact two-sided sign test (non-zero m={m_nz}): p = {p_sign:.6g}")

    # Exact within-(producer x week) reviewer-label permutation test
    by_pw: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_pw[(r["producer"], r["week"])].append(r)
    pw_overlap = sorted(k for k, rs in by_pw.items() if {r["reviewer"] for r in rs} == {"R1", "R2"})
    pw_tables: list[tuple[int, int, int, int]] = []
    for k in pw_overlap:
        rs = by_pw[k]
        r1s = [r for r in rs if r["reviewer"] == "R1"]
        r2s = [r for r in rs if r["reviewer"] == "R2"]
        a = sum(r["last_decision"] == "approve" for r in r1s)
        c = sum(r["last_decision"] == "approve" for r in r2s)
        pw_tables.append((a, len(r1s) - a, c, len(r2s) - c))
    pw_cmh = exact_cmh(pw_tables)
    print(
        f"  Exact within-(producer x week) permutation/CMH ({len(pw_overlap)} overlapping cells, "
        f"{sum(len(by_pw[k]) for k in pw_overlap)} rows): "
        f"S_obs={int(pw_cmh['s_obs'])} (support [{int(pw_cmh['s_min'])}, {int(pw_cmh['s_max'])}]) · "
        f"OR_MH={_fmt_or(pw_cmh['or_mh'])} · exact CMLE OR={_fmt_or(pw_cmh['cmle_or'])} "
        f"(95% CI [{_fmt_or(pw_cmh['exact_lo'])}, {_fmt_or(pw_cmh['exact_hi'])}]) · "
        f"exact two-sided p = {pw_cmh['p_two_sided']:.4e} · MH RD = {100.0 * pw_cmh['rd_mh']:+.1f} pp"
    )

    # Also Monte Carlo permutation of unweighted mean paired diff across (producer x week) cells
    rng_perm = random.Random(SEED)
    pw_cell_data: list[tuple[int, list[int]]] = []
    obs_cell_diffs: list[float] = []
    for k in pw_overlap:
        rs = by_pw[k]
        n1 = sum(r["reviewer"] == "R1" for r in rs)
        n2 = len(rs) - n1
        y = [1 if r["last_decision"] == "approve" else 0 for r in rs]
        a_obs = sum(1 for r in rs if r["reviewer"] == "R1" and r["last_decision"] == "approve")
        c_obs = sum(1 for r in rs if r["reviewer"] == "R2" and r["last_decision"] == "approve")
        obs_cell_diffs.append(a_obs / n1 - c_obs / n2)
        pw_cell_data.append((n1, y))
    obs_mean_pw = sum(obs_cell_diffs) / len(obs_cell_diffs)
    ge_abs = 0
    for _ in range(PERM_B):
        tot_d = 0.0
        for n1, y in pw_cell_data:
            n_tot = len(y)
            n2 = n_tot - n1
            # Sample n1 outcomes without replacement for R1
            r1_sum = sum(rng_perm.sample(y, n1))
            r2_sum = sum(y) - r1_sum
            tot_d += (r1_sum / n1) - (r2_sum / n2)
        if abs(tot_d / len(pw_cell_data)) >= abs(obs_mean_pw) - 1e-12:
            ge_abs += 1
    p_perm_mc = (ge_abs + 1) / (PERM_B + 1)
    print(
        f"  Monte Carlo within-(producer x week) permutation on unweighted mean cell diff "
        f"(B={PERM_B}, seed={SEED}): obs={100*obs_mean_pw:+.1f} pp, two-sided p = {p_perm_mc:.5f}"
    )

    # (c2) Producer-level cluster bootstrap 95% CI
    print(f"### (c2) 생산자 단위 Cluster Bootstrap 95% CI (B = {BOOT_B}, seed = {SEED})")

    def _cluster_boot(prod_ids: list[str], subset_rows: list[dict[str, str]], label: str, rng: random.Random) -> None:
        # Pre-aggregate per producer: (a1, n1, a2, n2)
        agg: dict[str, tuple[int, int, int, int]] = {}
        for p in prod_ids:
            rs = [r for r in subset_rows if r["producer"] == p]
            r1s = [r for r in rs if r["reviewer"] == "R1"]
            r2s = [r for r in rs if r["reviewer"] == "R2"]
            agg[p] = (
                sum(r["last_decision"] == "approve" for r in r1s),
                len(r1s),
                sum(r["last_decision"] == "approve" for r in r2s),
                len(r2s),
            )
        t_a1 = sum(v[0] for v in agg.values())
        t_n1 = sum(v[1] for v in agg.values())
        t_a2 = sum(v[2] for v in agg.values())
        t_n2 = sum(v[3] for v in agg.values())
        pt = (t_a1 / t_n1) - (t_a2 / t_n2)

        vals: list[float] = []
        k_prods = len(prod_ids)
        for _ in range(BOOT_B):
            sa1 = sn1 = sa2 = sn2 = 0
            for _ in range(k_prods):
                p = prod_ids[rng.randrange(k_prods)]
                a1, n1, a2, n2 = agg[p]
                sa1 += a1
                sn1 += n1
                sa2 += a2
                sn2 += n2
            if sn1 > 0 and sn2 > 0:
                vals.append((sa1 / sn1) - (sa2 / sn2))
        vals.sort()
        nv = len(vals)
        lo = vals[int(0.025 * nv)]
        hi = vals[min(nv - 1, int(0.975 * nv))]
        mean_v = sum(vals) / nv
        se_v = math.sqrt(sum((v - mean_v) ** 2 for v in vals) / (nv - 1))
        print(
            f"  {label:38} (clusters={k_prods:2}, valid={nv}/{BOOT_B}): "
            f"point={100*pt:+5.1f} pp · 95% CI [{100*lo:+5.1f}, {100*hi:+5.1f}] pp · SE={100*se_v:.1f} pp"
        )

    rng_boot = random.Random(SEED)
    all_prods = sorted(by_prod.keys())
    _cluster_boot(all_prods, rows, "All 21 producers (pooled R1-R2)", rng_boot)
    codex_prods = sorted({r["producer"] for r in rows if r["producer_tribe"] == "codex_cli"})
    _cluster_boot(codex_prods, [r for r in rows if r["producer_tribe"] == "codex_cli"], "codex_cli producers (all weeks)", rng_boot)
    agy_prods = sorted({r["producer"] for r in rows if r["producer_tribe"] == "agy_cli"})
    _cluster_boot(agy_prods, [r for r in rows if r["producer_tribe"] == "agy_cli"], "agy_cli producers (all weeks)", rng_boot)
    w4647_rows = [r for r in rows if r["week"] in ("2026-W46", "2026-W47") and r["producer_tribe"] in ("codex_cli", "agy_cli")]
    w4647_prods = sorted({r["producer"] for r in w4647_rows})
    _cluster_boot(w4647_prods, w4647_rows, "W46-W47 codex+agy overlap producers", rng_boot)
    _cluster_boot(both_prods, [r for r in rows if r["producer"] in set(both_prods)], "11 both-reviewer producers (pooled)", rng_boot)
    print()

    # ------------------------------------------------------------------
    # (d) Assignment imbalance — Review §3.2
    # ------------------------------------------------------------------
    print("## (d) 배정 쏠림 — W46·W47 생산자별 R1/R2 배정 동질성 검정 [Review §3.2]")
    print("  [boundary] 표에 있는 (week, producer_tribe, producer, reviewer) 배정 빈도만 검정함; 원인(크기 상한·과업 종류)은 주장하지 않음.")
    rng_assign = random.Random(SEED)
    for wk in ("2026-W46", "2026-W47"):
        for trb in ("codex_cli", "agy_cli"):
            sub = [r for r in rows if r["week"] == wk and r["producer_tribe"] == trb]
            prods = sorted({r["producer"] for r in sub})
            counts: list[tuple[int, int]] = []
            parts: list[str] = []
            for p in prods:
                r1c = sum(1 for r in sub if r["producer"] == p and r["reviewer"] == "R1")
                r2c = sum(1 for r in sub if r["producer"] == p and r["reviewer"] == "R2")
                counts.append((r1c, r2c))
                parts.append(f"{p}:{r1c}/{r1c+r2c}")
            chi2, p_val = perm_homogeneity_test(counts, n_perm=PERM_B, rng=rng_assign)
            tot_r1 = sum(x for x, _ in counts)
            tot_n = sum(x + y for x, y in counts)
            print(f"  {wk}  {trb:10}  total R1={tot_r1}/{tot_n} ({100*tot_r1/tot_n:.1f}%) · producers={len(prods)} · chi2={chi2:.3f} · exact-perm p={p_val:.5f}")
            print(f"    per-producer R1/total: {'  '.join(parts)}")
    # Also pooled W46+W47 by tribe
    for trb in ("codex_cli", "agy_cli"):
        sub = [r for r in rows if r["week"] in ("2026-W46", "2026-W47") and r["producer_tribe"] == trb]
        prods = sorted({r["producer"] for r in sub})
        counts = [
            (
                sum(1 for r in sub if r["producer"] == p and r["reviewer"] == "R1"),
                sum(1 for r in sub if r["producer"] == p and r["reviewer"] == "R2"),
            )
            for p in prods
        ]
        chi2, p_val = perm_homogeneity_test(counts, n_perm=PERM_B, rng=rng_assign)
        tot_r1 = sum(x for x, _ in counts)
        tot_n = sum(x + y for x, y in counts)
        parts = [f"{p}:{r1}/{r1+r2}" for p, (r1, r2) in zip(prods, counts)]
        print(f"  W46+W47  {trb:10}  total R1={tot_r1}/{tot_n} ({100*tot_r1/tot_n:.1f}%) · producers={len(prods)} · chi2={chi2:.3f} · exact-perm p={p_val:.5f}")
        print(f"    per-producer R1/total: {'  '.join(parts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
