#!/usr/bin/env python3
"""
order_distributions.py

For each chosen n, compute the multiplicative order of every unit in
(Z/nZ)* and plot a histogram of those orders. Peaks land at divisors of
lambda(n); the multiplicities are determined by the invariant factor
decomposition.

This is the parallel-cycle model made visible at the level of individual
elements. Cyclic groups produce one element per order divisor; fractured
groups produce many elements at the same low orders.

Generates: order_distributions.png
"""

from __future__ import annotations

from collections import Counter

import matplotlib.pyplot as plt

from lambda_ratio_explorer import (
    carmichael_lambda,
    collapse_index,
    element_orders,
    euler_totient,
    invariant_factors,
)


CASES: list[tuple[int, str]] = [
    (7,    "n = 7  (prime, fully cyclic)"),
    (24,   "n = 24 = 2^3 . 3  (maximally fractured small composite)"),
    (1729, "n = 1729 = 7 . 13 . 19  (Carmichael number)"),
    (1365, "n = 1365 = 3 . 5 . 7 . 13  (most-collapsed in [2, 2000])"),
]


def format_decomp(factors: list[int]) -> str:
    if factors == [1]:
        return "trivial"
    return " x ".join(f"Z/{d}" for d in factors)


def panel(ax: plt.Axes, n: int, title: str) -> None:
    orders = element_orders(n)
    counts = Counter(orders.values())

    phi = euler_totient(n)
    lam = carmichael_lambda(n)
    C = collapse_index(n)
    decomp = invariant_factors(n)

    keys = sorted(counts.keys())
    heights = [counts[k] for k in keys]

    bars = ax.bar(range(len(keys)), heights,
                  color="#4263eb", alpha=0.75, edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys)
    ax.set_xlabel("element order (must divide lambda(n))")
    ax.set_ylabel("count of units with that order")
    ax.set_title(title, fontsize=11)
    ax.grid(True, axis="y", alpha=0.2)

    if lam in keys:
        idx = keys.index(lam)
        bars[idx].set_color("#e8590c")

    info = (
        f"phi(n)    = {phi}\n"
        f"lambda(n) = {lam}\n"
        f"C(n)      = {C}\n"
        f"fracture  = {len(decomp)}\n"
        f"(Z/nZ)* = {format_decomp(decomp)}"
    )
    ax.text(0.97, 0.97, info,
            transform=ax.transAxes,
            fontsize=8.5,
            family="monospace",
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="white",
                      alpha=0.92, edgecolor="#adb5bd"))


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle(
        "Element-order distributions inside (Z/nZ)*  "
        "(orange bar = lambda(n), the longest stride)",
        fontsize=13, fontweight="bold",
    )

    for ax, (n, title) in zip(axes.flat, CASES):
        panel(ax, n, title)

    fig.tight_layout()
    out = "order_distributions.png"
    fig.savefig(out, dpi=160)
    print(f"Wrote {out}")

    print(f"\n{'-' * 60}")
    for n, title in CASES:
        orders = element_orders(n)
        counts = Counter(orders.values())
        decomp = invariant_factors(n)
        print(f"\n  n = {n}")
        print(f"    decomposition:  {format_decomp(decomp)}")
        print(f"    order -> count: "
              + ", ".join(f"{o}:{c}" for o, c in sorted(counts.items())))


if __name__ == "__main__":
    main()
