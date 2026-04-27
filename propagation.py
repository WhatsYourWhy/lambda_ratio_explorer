#!/usr/bin/env python3
"""
propagation.py

Demonstration of the multiplicative collapse propagation theorem

    C(k * p) = C(k) * gcd(lambda(k), p - 1)

for prime p coprime to k. The theorem is proved in THEOREM.md and
implemented as collapse_step / collapse_propagation_trace in
lambda_ratio_explorer.py.

This script:
  1. Builds three example numbers one prime at a time, printing the
     full propagation trace and showing how each new prime contributes
     its gcd factor to the running collapse index.
  2. Plots a 2-panel figure (propagation.png):
       - left:  step-by-step gcd contributions for 1365 and 1729
       - right: empirical density of gcd(p-1, q-1) over odd prime pairs
                up to a configurable limit (Option B-lite, no asymptotic
                claims; the figure simply shows where the mass sits).
"""

from __future__ import annotations

import math
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np

from lambda_ratio_explorer import collapse_propagation_trace


EXAMPLES: list[tuple[str, list[int]]] = [
    ("1365 = 3 * 5 * 7 * 13",                 [3, 5, 7, 13]),
    ("1729 = 7 * 13 * 19  (Carmichael)",      [7, 13, 19]),
    ("253  = 11 * 23   (clean shape)",        [11, 23]),
    ("91   = 7 * 13   (overlapping shape)",   [7, 13]),
]

GCD_HEATMAP_PRIME_LIMIT = 500


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i, b in enumerate(sieve) if b]


def print_trace(label: str, primes: list[int]) -> list[dict]:
    print(f"\n  {label}")
    trace = collapse_propagation_trace(primes)
    print(f"    {'step':>4}  {'k':>6}  {'p':>4}  {'lambda(k)':>10}  "
          f"{'gcd':>4}  {'-> kp':>8}  {'C(kp)':>6}")
    for i, step in enumerate(trace, start=1):
        print(f"    {i:4d}  {step['k']:6d}  {step['p']:4d}  "
              f"{step['lambda_k']:10d}  {step['gcd_term']:4d}  "
              f"{step['kp']:8d}  {step['C_kp']:6d}")
    return trace


def panel_propagation(ax: plt.Axes,
                      traces: list[tuple[str, list[dict]]]) -> None:
    """Stacked-bar visualization: each prime contributes a gcd factor.

    A vertical bar per example, broken into segments whose heights are
    log gcd_term. The total height equals log C(final).
    """
    labels = [t[0] for t in traces]
    n = len(traces)

    base_palette = ["#1971c2", "#ae3ec9", "#0ca678", "#fcc419", "#fd7e14",
                    "#e8590c", "#c92a2a"]

    bar_x = np.arange(n)
    bottom = np.zeros(n)

    max_steps = max(len(t[1]) for t in traces)
    for step_idx in range(max_steps):
        heights = []
        annotations = []
        for _, trace in traces:
            if step_idx < len(trace):
                step = trace[step_idx]
                heights.append(math.log2(step["gcd_term"]) if step["gcd_term"] > 1 else 0.0)
                annotations.append(f"x{step['p']}: gcd={step['gcd_term']}")
            else:
                heights.append(0.0)
                annotations.append("")
        heights_arr = np.array(heights)
        color = base_palette[step_idx % len(base_palette)]
        ax.bar(bar_x, heights_arr, bottom=bottom,
               color=color, edgecolor="white", linewidth=0.7,
               label=f"step {step_idx + 1}")
        for j, h in enumerate(heights_arr):
            if h > 0.05:
                ax.text(bar_x[j], bottom[j] + h / 2, annotations[j],
                        ha="center", va="center", fontsize=8,
                        color="white", fontweight="bold")
        bottom = bottom + heights_arr

    for j, (_, trace) in enumerate(traces):
        if not trace:
            continue
        final_C = trace[-1]["C_kp"]
        ax.text(bar_x[j], bottom[j] + 0.1, f"C = {final_C}",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(bar_x)
    ax.set_xticklabels(labels, fontsize=9, rotation=15, ha="right")
    ax.set_ylabel("log2 of multiplicative gcd contribution")
    ax.set_title("Step-by-step collapse propagation\n"
                 "each segment = log2 gcd(lambda(k), p - 1) for one new prime")
    ax.legend(fontsize="x-small", loc="upper left")
    ax.grid(True, axis="y", alpha=0.2)


def panel_gcd_distribution(ax: plt.Axes, prime_limit: int) -> None:
    """Empirical density of gcd(p-1, q-1) over distinct odd prime pairs."""
    primes = [p for p in primes_up_to(prime_limit) if p > 2]
    gcds: list[int] = []
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            gcds.append(math.gcd(p - 1, q - 1))

    counts = Counter(gcds)
    total = sum(counts.values())
    keys = sorted(counts.keys())
    densities = [counts[k] / total for k in keys]

    ax.bar(range(len(keys)), densities,
           color="#4263eb", alpha=0.8, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(keys, fontsize=8, rotation=0)
    ax.set_xlabel("gcd(p - 1, q - 1)")
    ax.set_ylabel("empirical density")
    ax.set_title(
        f"Density of gcd(p-1, q-1) over odd prime pairs, p,q <= {prime_limit}\n"
        f"({len(primes)} primes -> {total} pairs)"
    )
    ax.grid(True, axis="y", alpha=0.2)

    most_common = counts.most_common(3)
    summary = "  ".join(f"{k}: {v / total:.3f}" for k, v in most_common)
    ax.text(0.97, 0.95, "top values\n" + summary,
            transform=ax.transAxes,
            fontsize=8.5, family="monospace",
            verticalalignment="top", horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="white",
                      alpha=0.9, edgecolor="#adb5bd"))


def main() -> None:
    print("Collapse propagation traces")
    print("=" * 60)
    traces: list[tuple[str, list[dict]]] = []
    for label, primes in EXAMPLES:
        trace = print_trace(label, primes)
        traces.append((label, trace))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(
        "Collapse propagation:  C(k*p) = C(k) * gcd(lambda(k), p - 1)",
        fontsize=13, fontweight="bold",
    )
    panel_propagation(axes[0], traces)
    panel_gcd_distribution(axes[1], GCD_HEATMAP_PRIME_LIMIT)

    fig.tight_layout()
    out = "propagation.png"
    fig.savefig(out, dpi=160)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
