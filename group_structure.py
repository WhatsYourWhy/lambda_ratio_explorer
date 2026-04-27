#!/usr/bin/env python3
"""
group_structure.py

Visualize multiplicative-group collapse via the index

    C(n) = phi(n) / lambda(n)

For semiprimes n = p*q this reduces to gcd(p-1, q-1), which is exactly the
quantity controlling RSA's resistance to small-subgroup and short-cycle
attacks. For general n, C(n) counts how many cyclic factors the group
(Z/nZ)* decomposes into beyond its largest one.

Generates a 4-panel teaching figure: group_structure.png
  1. C(n) across q   -- structural collapse vs. n
  2. Semiprime heatmap of log2 gcd(p-1, q-1)  -- the cryptographic map
  3. Distribution of C(n) by structural kind
  4. lambda(n)/(n-1) vs. 1/C(n)  -- two views of the same collapse
"""

from __future__ import annotations

import math
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from lambda_ratio_explorer import (
    carmichael_lambda,
    collapse_index,
    euler_totient,
    factorize,
    is_carmichael,
    is_prime,
    kind,
)

Q_MAX = 2000
HEATMAP_PRIMES_UP_TO = 400


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i, b in enumerate(sieve) if b]


def panel_collapse_scan(ax: plt.Axes, q_max: int) -> None:
    """Panel 1: C(q) vs q, colored by structural kind."""
    qs = list(range(2, q_max + 1))
    kinds = [kind(q) for q in qs]
    cs = [collapse_index(q) for q in qs]

    by_kind: dict[str, tuple[list[int], list[int]]] = {
        "prime": ([], []),
        "prime_power": ([], []),
        "composite": ([], []),
    }
    for q, k, c in zip(qs, kinds, cs):
        by_kind[k][0].append(q)
        by_kind[k][1].append(c)

    style = {
        "prime":       {"color": "#e8590c", "s": 10, "marker": "o", "alpha": 0.6, "zorder": 4},
        "prime_power": {"color": "#2b8a3e", "s": 22, "marker": "s", "alpha": 0.6, "zorder": 3},
        "composite":   {"color": "#4263eb", "s":  6, "marker": ".", "alpha": 0.4, "zorder": 2},
    }
    for k in ("composite", "prime_power", "prime"):
        xs, ys = by_kind[k]
        ax.scatter(xs, ys, label=k, **style[k])

    carmichaels = [(q, c) for q, c in zip(qs, cs) if is_carmichael(q)]
    if carmichaels:
        cx, cy = zip(*carmichaels)
        ax.scatter(cx, cy, s=70, marker="*", color="#fab005",
                   edgecolor="black", linewidth=0.6, zorder=5,
                   label=f"Carmichael number ({len(carmichaels)} found)")

    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("q")
    ax.set_ylabel("C(q) = phi(q) / lambda(q)")
    ax.set_title("Collapse index across the integers")
    ax.axhline(1, color="black", ls="--", alpha=0.3, linewidth=0.8,
               label="C = 1 (cyclic group)")
    ax.legend(fontsize="small", loc="upper left")
    ax.grid(True, alpha=0.2, which="both")


def panel_semiprime_heatmap(ax: plt.Axes, prime_limit: int) -> None:
    """Panel 2: heatmap of log2 gcd(p-1, q-1) over odd prime pairs."""
    primes = [p for p in primes_up_to(prime_limit) if p > 2]
    n = len(primes)
    grid = np.zeros((n, n), dtype=float)

    for i, p in enumerate(primes):
        for j, q in enumerate(primes):
            if i == j:
                grid[i, j] = np.nan
                continue
            grid[i, j] = math.gcd(p - 1, q - 1)

    valid = grid[~np.isnan(grid)]
    vmin = max(1.0, valid.min())
    vmax = valid.max()

    im = ax.imshow(
        grid,
        origin="lower",
        cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        extent=(primes[0], primes[-1], primes[0], primes[-1]),
        aspect="equal",
        interpolation="nearest",
    )
    ax.set_xlabel("p")
    ax.set_ylabel("q")
    ax.set_title(
        f"C(pq) = gcd(p-1, q-1)  for odd primes up to {prime_limit}\n"
        "dark = strong group structure (good for RSA)"
    )
    cbar = plt.colorbar(im, ax=ax, label="gcd(p-1, q-1)", shrink=0.85)
    cbar.ax.tick_params(labelsize="x-small")


def panel_distribution(ax: plt.Axes, q_max: int) -> None:
    """Panel 3: distribution of C(q) by kind."""
    by_kind: dict[str, list[int]] = {"prime": [], "prime_power": [], "composite": []}
    for q in range(2, q_max + 1):
        by_kind[kind(q)].append(collapse_index(q))

    max_C = max(max(v) for v in by_kind.values() if v)
    bins = np.unique(np.geomspace(1, max_C + 1, num=40).round().astype(int))

    colors = {"prime": "#e8590c", "prime_power": "#2b8a3e", "composite": "#4263eb"}
    for k in ("prime", "prime_power", "composite"):
        vals = by_kind[k]
        if not vals:
            continue
        ax.hist(vals, bins=bins, alpha=0.55, color=colors[k],
                label=f"{k} (n={len(vals)})", edgecolor="white", linewidth=0.4)

    ax.set_xscale("log")
    ax.set_xlabel("C(q)")
    ax.set_ylabel("count")
    ax.set_title("Distribution of collapse indices")
    ax.legend(fontsize="small")
    ax.grid(True, alpha=0.2, which="both")


def panel_two_views(ax: plt.Axes, q_max: int) -> None:
    """Panel 4: lambda(q)/(q-1) vs. 1/C(q).

    Both quantities measure how cyclic the group (Z/qZ)* is.
    Plotting them together shows where they agree and where they diverge.
    """
    qs, lam_ratio, inv_c, kinds = [], [], [], []
    for q in range(3, q_max + 1):
        lam = carmichael_lambda(q)
        c = collapse_index(q)
        qs.append(q)
        lam_ratio.append(lam / (q - 1))
        inv_c.append(1 / c)
        kinds.append(kind(q))

    style = {
        "prime":       {"color": "#e8590c", "s": 10, "marker": "o", "alpha": 0.7},
        "prime_power": {"color": "#2b8a3e", "s": 22, "marker": "s", "alpha": 0.7},
        "composite":   {"color": "#4263eb", "s":  6, "marker": ".", "alpha": 0.35},
    }

    for k in ("composite", "prime_power", "prime"):
        idx = [i for i, ki in enumerate(kinds) if ki == k]
        ax.scatter(
            [lam_ratio[i] for i in idx],
            [inv_c[i] for i in idx],
            label=k, **style[k],
        )

    diag = np.linspace(0, 1.05, 50)
    ax.plot(diag, diag, "k--", alpha=0.4, linewidth=0.8, label="y = x")

    ax.set_xlabel("lambda(q) / (q - 1)")
    ax.set_ylabel("1 / C(q)   [cyclic fraction]")
    ax.set_title("Two collapse measures, side by side")
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize="small", loc="lower right")
    ax.grid(True, alpha=0.2)


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle(
        "Multiplicative-group collapse:  C(n) = phi(n) / lambda(n)",
        fontsize=14, fontweight="bold",
    )

    panel_collapse_scan(axes[0, 0], Q_MAX)
    panel_semiprime_heatmap(axes[0, 1], HEATMAP_PRIMES_UP_TO)
    panel_distribution(axes[1, 0], Q_MAX)
    panel_two_views(axes[1, 1], Q_MAX)

    fig.tight_layout()
    out = "group_structure.png"
    fig.savefig(out, dpi=160)
    print(f"Wrote {out}")

    print_summary(Q_MAX)


def print_summary(q_max: int) -> None:
    by_kind: dict[str, list[int]] = {"prime": [], "prime_power": [], "composite": []}
    for q in range(2, q_max + 1):
        by_kind[kind(q)].append(collapse_index(q))

    print(f"\n{'-' * 60}")
    print(f"  C(q) = phi(q) / lambda(q)   for q in [2, {q_max}]")
    print(f"{'-' * 60}")
    for k in ("prime", "prime_power", "composite"):
        vals = np.array(by_kind[k])
        if len(vals) == 0:
            continue
        print(f"\n  {k} (count={len(vals)}):")
        print(f"    mean C   = {vals.mean():.3f}")
        print(f"    median C = {int(np.median(vals))}")
        print(f"    max C    = {int(vals.max())}")
        print(f"    fraction with C=1: {(vals == 1).mean():.3f}")

    extremes = sorted(
        ((q, collapse_index(q)) for q in range(2, q_max + 1)),
        key=lambda x: x[1], reverse=True,
    )[:10]
    print("\n  Top 10 most-collapsed q (smallest cyclic fraction):")
    for q, c in extremes:
        f = factorize(q)
        fstr = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))
        print(f"    q={q:5d}  C={c:5d}   factors: {fstr}")


if __name__ == "__main__":
    main()
