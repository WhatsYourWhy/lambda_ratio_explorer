#!/usr/bin/env python3
"""
gcd_distribution_theory.py

Empirical verification of Theorem C from THEOREM.md.

Theorem C: for distinct odd primes p, q drawn uniformly from primes <= X,
and any prime l,

    Pr( l | gcd(p - 1, q - 1) )  ->  1 / (l - 1)^2  as X -> infty.

This script:
  1. Sieves odd primes up to a configurable limit.
  2. For l in {2, 3, 5, 7, 11, 13, 17, 19, 23} computes empirical
     Pr( l | gcd(p - 1, q - 1) ) over all distinct odd prime pairs and
     compares to the Dirichlet prediction 1/(l-1)^2.
  3. Computes empirical E[gcd(p - 1, q - 1)] and the truncated heuristic
     sum_{d <= D} 1/phi(d), which grows like A log D with
     A = 315 zeta(3) / (2 pi^4).
  4. Renders a two-panel figure gcd_distribution.png:
       - left:  empirical vs predicted Pr(l | gcd) on a log y-axis
       - right: the 1365 = 3*5*7*13 tie-back, showing how its propagation
                trace concentrates exactly the small-prime gcd factors
                that Theorem C predicts dominate the distribution.

The script prints residuals so the strength of the agreement is visible
at a glance.
"""

from __future__ import annotations

import math
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from lambda_ratio_explorer import (
    collapse_propagation_trace,
    euler_totient,
)


PRIME_LIMIT = 2000
TARGET_PRIMES: list[int] = [2, 3, 5, 7, 11, 13, 17, 19, 23]

# Mean-totient-reciprocal constant A = 315 zeta(3) / (2 pi^4),
# governing the asymptotic sum_{d <= X} 1/phi(d) ~ A log X.
ASYMPTOTIC_CONSTANT_A = 1.9436


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


def empirical_divisibility_rates(
    primes: list[int], target_ells: Iterable[int]
) -> tuple[dict[int, float], int, float]:
    """Return (empirical Pr(l | gcd) per l, total pair count, empirical mean gcd)."""
    counts = {l: 0 for l in target_ells}
    pair_total = 0
    gcd_sum = 0
    for i, p in enumerate(primes):
        for q in primes[i + 1:]:
            g = math.gcd(p - 1, q - 1)
            pair_total += 1
            gcd_sum += g
            for l in counts:
                if g % l == 0:
                    counts[l] += 1
    rates = {l: counts[l] / pair_total for l in counts}
    mean_gcd = gcd_sum / pair_total
    return rates, pair_total, mean_gcd


def dirichlet_predicted_rate(l: int) -> float:
    return 1.0 / (l - 1) ** 2


def truncated_mean_heuristic(D: int) -> float:
    """sum_{d <= D} 1/phi(d), the Dirichlet-rate substitution into the
    standard identity gcd(a, b) = sum_d phi(d) [d|a][d|b]."""
    return sum(1.0 / euler_totient(d) for d in range(1, D + 1))


def print_table(
    rates: dict[int, float],
    pair_total: int,
    mean_gcd: float,
    truncated_heuristic: float,
    asymptotic_estimate: float,
    prime_limit: int,
) -> None:
    print("\nTheorem C: empirical vs Dirichlet-predicted Pr(l | gcd(p-1, q-1))")
    print(f"  over {pair_total} distinct odd-prime pairs")
    print()
    print(f"   {'l':>4}  {'empirical':>10}  {'predicted':>10}  {'abs error':>10}  {'rel error':>10}")
    print(f"   {'-'*4:>4}  {'-'*10:>10}  {'-'*10:>10}  {'-'*10:>10}  {'-'*10:>10}")
    for l, emp in rates.items():
        pred = dirichlet_predicted_rate(l)
        abs_err = emp - pred
        rel_err = abs_err / pred if pred > 0 else 0.0
        print(f"   {l:4d}  {emp:10.4f}  {pred:10.4f}  {abs_err:+10.4f}  {rel_err:+10.2%}")

    print()
    print("Expected gcd, asymptotic E[gcd] ~ A log X with A = 315 zeta(3) / (2 pi^4) ~= 1.9436:")
    print(f"   empirical mean over pairs                           = {mean_gcd:8.4f}")
    print(f"   sum_{{d <= {prime_limit}}} 1 / phi(d)  (Dirichlet substitution) = {truncated_heuristic:8.4f}")
    print(f"   A * log({prime_limit})                                       = {asymptotic_estimate:8.4f}")
    print()
    print("Empirical falls below the heuristic at finite X because Pr(d | p-1)")
    print("does not yet match 1/phi(d) for d comparable to X (slow convergence")
    print("of Dirichlet density for large moduli). The takeaway: E[gcd] grows")
    print("logarithmically with the prime size, so safe-prime selection (which")
    print("caps the gcd at 2 regardless of X) is the practical RSA defense.")


def panel_rates(
    ax: plt.Axes,
    rates: dict[int, float],
    target_ells: list[int],
    pair_total: int,
) -> None:
    ells = list(target_ells)
    n = len(ells)
    width = 0.38
    x = np.arange(n)

    empirical = [rates[l] for l in ells]
    predicted = [dirichlet_predicted_rate(l) for l in ells]

    ax.bar(x - width / 2, empirical, width=width,
           label="empirical", color="#4263eb",
           edgecolor="white", linewidth=0.5)
    ax.bar(x + width / 2, predicted, width=width,
           label="$1/(\\ell-1)^2$ (Dirichlet)", color="#fa5252",
           alpha=0.85, edgecolor="white", linewidth=0.5)

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([str(l) for l in ells])
    ax.set_xlabel("prime $\\ell$")
    ax.set_ylabel("Pr$(\\ell \\mid \\gcd(p - 1, q - 1))$  (log scale)")
    ax.set_title(
        f"Theorem C: divisibility of $\\gcd(p-1, q-1)$\n"
        f"odd-prime pairs from primes $\\leq$ {PRIME_LIMIT}, "
        f"{pair_total} pairs"
    )
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(loc="upper right", fontsize="small")

    ax.annotate(
        "$\\ell = 2$ saturates: every odd\nprime has $p - 1$ even",
        xy=(0, 1.0), xytext=(0.6, 0.35),
        fontsize=8.5,
        arrowprops=dict(arrowstyle="->", color="#495057", lw=0.8),
        color="#495057",
    )


def panel_1365_tieback(ax: plt.Axes) -> None:
    """The 1365 = 3*5*7*13 propagation, colored by prime contribution.

    Theorem C says small primes carry the mass of the gcd distribution.
    1365's prime factors 3, 5, 7, 13 produce totients 2, 4, 6, 12, all
    rich in small-prime structure. The propagation trace shows each step
    cashing in exactly the small-prime overlap Theorem C predicts.
    """
    primes = [3, 5, 7, 13]
    trace = collapse_propagation_trace(primes)

    palette = ["#1971c2", "#ae3ec9", "#0ca678", "#fcc419"]

    bottom = 0.0
    bar_x = 0
    for step_idx, step in enumerate(trace):
        gcd_term = step["gcd_term"]
        height = math.log2(gcd_term) if gcd_term > 1 else 0.0
        color = palette[step_idx % len(palette)]
        ax.bar(bar_x, height, bottom=bottom, width=0.55,
               color=color, edgecolor="white", linewidth=0.7,
               label=f"x{step['p']}: gcd={gcd_term}")
        if height > 0.05:
            ax.text(bar_x, bottom + height / 2,
                    f"x{step['p']}\ngcd={gcd_term}",
                    ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold")
        bottom += height

    final_C = trace[-1]["C_kp"]
    ax.text(bar_x, bottom + 0.15, f"$C(1365) = {final_C}$",
            ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_xticks([bar_x])
    ax.set_xticklabels(["$1365 = 3 \\cdot 5 \\cdot 7 \\cdot 13$"])
    ax.set_xlim(-1, 1.4)
    ax.set_ylabel("$\\log_2$ gcd contribution per step")
    ax.set_title(
        "Why $1365$ is a Theorem C extremum\n"
        "small primes' totients $\\{2, 4, 6, 12\\}$ are exactly what\n"
        "Theorem C predicts dominates $\\gcd(p - 1, q - 1)$"
    )
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(loc="center right", fontsize="x-small",
              bbox_to_anchor=(1.0, 0.5))

    annotation = (
        "Theorem C predicts:\n"
        "  Pr$(2 \\mid \\gcd) = 1$\n"
        "  Pr$(3 \\mid \\gcd) = 1/4$\n"
        "  Pr$(4 \\mid \\gcd) = 1/4$\n"
        "1365's primes hit all three."
    )
    ax.text(0.02, 0.98, annotation, transform=ax.transAxes,
            fontsize=8.5, family="monospace",
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white",
                      alpha=0.9, edgecolor="#adb5bd"))


def main() -> None:
    primes = [p for p in primes_up_to(PRIME_LIMIT) if p > 2]
    rates, pair_total, mean_gcd = empirical_divisibility_rates(primes, TARGET_PRIMES)
    truncated_heuristic = truncated_mean_heuristic(PRIME_LIMIT)
    asymptotic_estimate = ASYMPTOTIC_CONSTANT_A * math.log(PRIME_LIMIT)
    print_table(
        rates, pair_total, mean_gcd, truncated_heuristic, asymptotic_estimate, PRIME_LIMIT
    )

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(
        "Theorem C: Dirichlet density of $\\gcd(p - 1, q - 1)$",
        fontsize=13, fontweight="bold",
    )
    panel_rates(axes[0], rates, TARGET_PRIMES, pair_total)
    panel_1365_tieback(axes[1])

    fig.tight_layout()
    out = "gcd_distribution.png"
    fig.savefig(out, dpi=160)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
