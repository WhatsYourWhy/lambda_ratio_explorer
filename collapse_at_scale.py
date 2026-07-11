#!/usr/bin/env python3
"""
collapse_at_scale.py

The collapse index C(n) = phi(n)/lambda(n) scanned to n = 10^6.

The per-n factorization in lambda_ratio_explorer.py is the right tool for
individual (possibly huge) n, but a full-range scan is faster with a
smallest-prime-factor sieve: one pass builds the factorization of every
n <= N, from which phi and lambda follow directly.

Generates a 2-panel figure: collapse_at_scale.png
  1. Density of (n, C(n)) over all composites on log-log axes, with the
     running-maximum "collapse records" overlaid and labeled. The record
     holders are exactly the products of small primes with heavily shared
     totient structure that Theorem C predicts.
  2. Quantile fan of C(n) by scale: median, 90th, 99th percentile, and max
     of C over geometric windows in n. The median stays in single digits
     while the extremes explode -- collapse is a tail phenomenon.

Also prints the by-kind summary table (the 10^6 version of the README
table) and the top collapse records with factorizations.
"""

from __future__ import annotations

import math
from array import array

import matplotlib.pyplot as plt
import numpy as np

from lambda_ratio_explorer import factorize

N_MAX = 1_000_000
N_QUANTILE_BINS = 40

# Repo palette: blues carry magnitude (sequential, light -> dark),
# orange marks the record trace, grays are ink.
BLUE_SEQ = ["#74c0fc", "#339af0", "#1971c2", "#0b4a8c"]
RECORD_COLOR = "#e8590c"
INK = "#495057"


def spf_sieve(limit: int) -> array:
    """Smallest prime factor of every n <= limit."""
    spf = array("l", range(limit + 1))
    for i in range(2, int(limit ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def scan(limit: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (C, kind_code) arrays indexed by n.

    kind_code: 0 = unit/unused, 1 = prime, 2 = prime power, 3 = composite.
    """
    spf = spf_sieve(limit)
    C = np.zeros(limit + 1, dtype=np.int64)
    kind_code = np.zeros(limit + 1, dtype=np.int8)
    C[1] = 1

    for n in range(2, limit + 1):
        m = n
        phi = 1
        lam = 1
        distinct = 0
        while m > 1:
            p = spf[m]
            k = 0
            while m % p == 0:
                m //= p
                k += 1
            distinct += 1
            phi *= (p - 1) * p ** (k - 1)
            if p == 2 and k >= 3:
                lam_pk = 2 ** (k - 2)
            else:
                lam_pk = (p - 1) * p ** (k - 1)
            lam = lam * lam_pk // math.gcd(lam, lam_pk)
            if distinct == 1 and m == 1:
                kind_code[n] = 1 if k == 1 else 2
        if kind_code[n] == 0:
            kind_code[n] = 3
        C[n] = phi // lam
    return C, kind_code


def collapse_records(C: np.ndarray) -> list[tuple[int, int]]:
    """(n, C(n)) points where C sets a new running maximum."""
    records = []
    best = 0
    for n in range(2, len(C)):
        if C[n] > best:
            best = int(C[n])
            records.append((n, best))
    return records


def print_summary(C: np.ndarray, kind_code: np.ndarray) -> None:
    names = {1: "prime", 2: "prime_power", 3: "composite"}
    print(f"\nCollapse index C(n) = phi(n)/lambda(n) over n in [2, {len(C) - 1}]\n")
    print(f"  {'kind':<12} {'count':>8} {'mean C':>9} {'median C':>9} {'max C':>8} {'frac C=1':>9}")
    for code, name in names.items():
        vals = C[2:][kind_code[2:] == code]
        print(f"  {name:<12} {len(vals):>8} {vals.mean():>9.2f} {int(np.median(vals)):>9} "
              f"{vals.max():>8} {(vals == 1).mean():>9.3f}")

    print("\nTop collapse records:")
    records = collapse_records(C)
    for n, c in records[-10:]:
        factors = " * ".join(
            (f"{p}^{k}" if k > 1 else str(p)) for p, k in factorize(n).items()
        )
        print(f"  C({n:>7}) = {c:>5}   {n} = {factors}")


def panel_density(ax: plt.Axes, C: np.ndarray, kind_code: np.ndarray) -> None:
    comp = np.flatnonzero(kind_code == 3)
    x = np.log10(comp.astype(np.float64))
    y = np.log2(C[comp].astype(np.float64))

    hb = ax.hexbin(x, y, gridsize=60, cmap="Blues", bins="log", mincnt=1,
                   linewidths=0.1)
    cb = plt.colorbar(hb, ax=ax, pad=0.01)
    cb.set_label("composites per cell (log scale)", fontsize=8.5)

    records = collapse_records(C)
    rx = [math.log10(n) for n, _ in records]
    ry = [math.log2(c) for _, c in records]
    ax.plot(rx, ry, color=RECORD_COLOR, lw=1.4, zorder=3,
            drawstyle="steps-post", label="running max of $C$")
    ax.scatter(rx, ry, color=RECORD_COLOR, s=22, zorder=4,
               edgecolor="white", linewidth=0.6)
    # Label a spread of records: the final one plus two mid-scale ones,
    # so annotations never pile up in the top-right corner.
    labeled = {records[-1]}
    for frac in (0.5, 0.78):
        target = frac * math.log10(records[-1][0])
        labeled.add(min(records, key=lambda r: abs(math.log10(r[0]) - target)))
    for n, c in labeled:
        ax.annotate(f"$C({n}) = {c}$", (math.log10(n), math.log2(c)),
                    xytext=(-8, 6), textcoords="offset points",
                    ha="right", fontsize=8.5, color=INK)

    xticks = range(1, 7)
    ax.set_xticks(list(xticks))
    ax.set_xticklabels([f"$10^{k}$" for k in xticks])
    yticks = range(0, int(max(ry)) + 2, 2)
    ax.set_yticks(list(yticks))
    ax.set_yticklabels([f"{2 ** k}" for k in yticks])
    ax.set_xlabel("$n$")
    ax.set_ylabel("$C(n)$")
    ax.set_title("Collapse of composites to $n = 10^6$\n"
                 "density of $(n, C)$ with record collapses overlaid")
    ax.grid(True, alpha=0.15)
    ax.legend(loc="upper left", fontsize="small")


def panel_quantile_fan(ax: plt.Axes, C: np.ndarray, kind_code: np.ndarray) -> None:
    comp = np.flatnonzero(kind_code == 3)
    Cc = C[comp].astype(np.float64)
    edges = np.geomspace(comp[0], len(C) - 1, N_QUANTILE_BINS + 1)

    quantiles = [(0.50, "median", BLUE_SEQ[0]),
                 (0.90, "90th pct", BLUE_SEQ[1]),
                 (0.99, "99th pct", BLUE_SEQ[2]),
                 (1.00, "max", BLUE_SEQ[3])]
    centers = np.sqrt(edges[:-1] * edges[1:])
    for q, label, color in quantiles:
        ys = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            window = Cc[(comp >= lo) & (comp < hi)]
            ys.append(np.quantile(window, q) if len(window) else np.nan)
        ax.plot(centers, ys, color=color, lw=1.8)
        ax.annotate(label, (centers[-1], ys[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center",
                    fontsize=8.5, color=INK)

    ax.set_xscale("log")
    ax.set_yscale("log", base=2)
    ax.set_xlabel("$n$ (geometric windows)")
    ax.set_ylabel("$C(n)$ quantile within window")
    ax.set_title("Collapse is a tail phenomenon\n"
                 "the median composite barely collapses; the extremes explode")
    ax.grid(True, which="both", alpha=0.15)
    ax.set_xlim(centers[0], centers[-1] * 3.2)


def main() -> None:
    print(f"sieving and scanning to {N_MAX} ...")
    C, kind_code = scan(N_MAX)
    print_summary(C, kind_code)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6.2))
    fig.suptitle("Collapse index $C(n) = \\varphi(n)/\\lambda(n)$ at scale",
                 fontsize=13, fontweight="bold")
    panel_density(axes[0], C, kind_code)
    panel_quantile_fan(axes[1], C, kind_code)

    fig.tight_layout()
    out = "collapse_at_scale.png"
    fig.savefig(out, dpi=160)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
