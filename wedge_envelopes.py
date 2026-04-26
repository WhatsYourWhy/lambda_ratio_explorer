#!/usr/bin/env python3
"""
wedge_envelopes.py

Companion visualization showing that the upper edges of the lambda(q) wedge
are exact algebraic lines, not approximations.

For semiprimes q = k*p with smallest factor k:
    lambda(k*p) = lcm(lambda(k), p - 1)
For k = 2:  lambda(2p) = p - 1            -> envelope (q/2 - 1)
For k = 3:  lambda(3p) = lcm(2, p-1) = p-1 (p > 2)  -> envelope (q/3 - 1)

Plots lambda(q) vs q with the predicted envelopes overlaid.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from lambda_ratio_explorer import carmichael_lambda, factorize, is_prime, kind

Q_MAX = 2000


def smallest_factor(q: int) -> int:
    return min(factorize(q))


def main() -> None:
    qs = list(range(2, Q_MAX + 1))
    lams = [carmichael_lambda(q) for q in qs]
    kinds = [kind(q) for q in qs]
    smallest = [smallest_factor(q) for q in qs]

    fig, ax = plt.subplots(figsize=(11, 7))

    prime_idx = [i for i, k in enumerate(kinds) if k == "prime"]
    ax.scatter([qs[i] for i in prime_idx], [lams[i] for i in prime_idx],
               s=10, color="#e8590c", alpha=0.6, label="prime")

    palette = {2: "#1971c2", 3: "#ae3ec9", 5: "#0ca678", 7: "#fcc419"}
    for k, color in palette.items():
        idx = [i for i, (ki, sf) in enumerate(zip(kinds, smallest))
               if ki == "composite" and sf == k]
        ax.scatter([qs[i] for i in idx], [lams[i] for i in idx],
                   s=12, color=color, alpha=0.55, label=f"smallest factor = {k}")

    other_idx = [i for i, (ki, sf) in enumerate(zip(kinds, smallest))
                 if ki == "composite" and sf not in palette]
    ax.scatter([qs[i] for i in other_idx], [lams[i] for i in other_idx],
               s=4, color="#ced4da", alpha=0.4, label="composite (other)")

    q_line = np.linspace(2, Q_MAX, 600)
    ax.plot(q_line, q_line - 1, "--", color="#e8590c", alpha=0.7,
            linewidth=1.4, label="prime ceiling: q - 1")
    for k, color in palette.items():
        ax.plot(q_line, q_line / k - 1, "--", color=color, alpha=0.7,
                linewidth=1.2, label=f"k=p envelope: q/{k} - 1")

    ax.set_xlabel("q")
    ax.set_ylabel("lambda(q)")
    ax.set_title("lambda(q) wedge: exact envelopes by smallest prime factor")
    ax.legend(fontsize="small", ncols=2, loc="upper left")
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    out = "wedge_envelopes.png"
    fig.savefig(out, dpi=160)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
