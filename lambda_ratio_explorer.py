#!/usr/bin/env python3
"""
lambda_ratio_explorer.py

Core library for studying multiplicative-group structure of (Z/nZ)*.

Computes three primary quantities for each n:

    phi(n)              -- Euler totient: order of (Z/nZ)*
    lambda(n)           -- Carmichael function: exponent of (Z/nZ)*
    C(n) = phi/lambda   -- collapse index: index of the cyclic factor

Interpretation:
    C(n) = 1            <=>  group is cyclic (n is prime, 2p, p^k, etc.)
    C(n) > 1            <=>  group decomposes; cycles overlap; structure collapses
    For n = p*q:        C(n) = gcd(p-1, q-1)
                        -- exactly the cryptographic weakness measure for RSA.

The optional ratio R(q, n) = lambda(q) / log(n) is retained for backward
compatibility but is not the central object: log(n) is just a scalar.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from functools import reduce
from math import gcd, lcm
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Row:
    q: int
    n: int
    kind: str
    phi_q: int
    lambda_q: int
    collapse: int
    ratio: float


def carmichael_lambda_bruteforce(n: int) -> int:
    """Small/slow definition check: smallest m where a^m = 1 mod n for all gcd(a,n)=1."""
    if n < 1:
        raise ValueError("n must be positive")
    if n == 1:
        return 1

    coprimes = [a for a in range(1, n) if gcd(a, n) == 1]
    for m in range(1, n * n + 1):
        if all(pow(a, m, n) == 1 for a in coprimes):
            return m
    raise RuntimeError(f"bruteforce search failed for n={n}")


def factorize(n: int) -> dict[int, int]:
    """Trial-division factorization, good enough for toy scans up to moderate q."""
    if n < 1:
        raise ValueError("n must be positive")

    factors: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2  # 2, then odd candidates only
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def carmichael_prime_power(p: int, k: int) -> int:
    """Carmichael lambda for p^k."""
    if p == 2 and k >= 3:
        return 2 ** (k - 2)
    return (p - 1) * (p ** (k - 1))


def carmichael_lambda(n: int, *, method: str = "fast") -> int:
    if method == "brute":
        return carmichael_lambda_bruteforce(n)
    if n == 1:
        return 1
    parts = [carmichael_prime_power(p, k) for p, k in factorize(n).items()]
    return reduce(lcm, parts, 1)


def euler_totient(n: int) -> int:
    """Euler's totient phi(n) = order of the multiplicative group (Z/nZ)*."""
    if n < 1:
        raise ValueError("n must be positive")
    if n == 1:
        return 1
    result = 1
    for p, k in factorize(n).items():
        result *= p ** (k - 1) * (p - 1)
    return result


def collapse_index(n: int) -> int:
    """C(n) = phi(n) / lambda(n).

    Equals the index of the largest cyclic factor in the decomposition of
    (Z/nZ)*. C(n) = 1 iff the group is cyclic. For n = p*q with p, q distinct
    odd primes, C(n) = gcd(p-1, q-1).
    """
    if n < 1:
        raise ValueError("n must be positive")
    if n == 1:
        return 1
    return euler_totient(n) // carmichael_lambda(n)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def kind(q: int) -> str:
    if is_prime(q):
        return "prime"
    if len(factorize(q)) == 1:
        return "prime_power"
    return "composite"


def ratio(q: int, n: int, *, method: str = "fast") -> float:
    if n <= 1:
        raise ValueError("n must be > 1 because log(n) is the normalizer")
    return carmichael_lambda(q, method=method) / math.log(n)


def scan(q_min: int, q_max: int, n_values: Iterable[int], *, method: str = "fast") -> list[Row]:
    rows: list[Row] = []
    for n in n_values:
        if n <= 1:
            raise ValueError("all n values must be > 1")
        log_n = math.log(n)
        for q in range(q_min, q_max + 1):
            lam = carmichael_lambda(q, method=method)
            phi = euler_totient(q)
            rows.append(Row(
                q=q, n=n, kind=kind(q),
                phi_q=phi, lambda_q=lam,
                collapse=phi // lam,
                ratio=lam / log_n,
            ))
    return rows


def print_table(rows: list[Row], *, limit: int | None = None) -> None:
    header = f"{'q':>5} {'n':>8} {'kind':>12} {'phi(q)':>8} {'lambda(q)':>10} {'C(q)':>6} {'lambda/log(n)':>15}"
    print(header)
    print("-" * len(header))
    for row in rows[:limit]:
        print(f"{row.q:5d} {row.n:8d} {row.kind:>12} {row.phi_q:8d} {row.lambda_q:10d} "
              f"{row.collapse:6d} {row.ratio:15.8f}")


def write_csv(rows: list[Row], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["q", "n", "kind", "phi_q", "lambda_q", "collapse", "ratio"]
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def summarize(rows: list[Row], *, top: int = 10) -> None:
    print("\nTop collapse indices C(q) = phi(q)/lambda(q):")
    for row in sorted(rows, key=lambda r: r.collapse, reverse=True)[:top]:
        print(f"  q={row.q:<6} {row.kind:<12} phi={row.phi_q:<6} lambda={row.lambda_q:<6} "
              f"C={row.collapse}")

    by_kind: dict[str, list[Row]] = {}
    for row in rows:
        by_kind.setdefault(row.kind, []).append(row)

    print("\nAverages by kind:")
    for k, items in sorted(by_kind.items()):
        avg_phi = sum(r.phi_q for r in items) / len(items)
        avg_lam = sum(r.lambda_q for r in items) / len(items)
        avg_C = sum(r.collapse for r in items) / len(items)
        print(f"  {k:<12} count={len(items):<5} avg_phi={avg_phi:.2f} "
              f"avg_lambda={avg_lam:.2f} avg_C={avg_C:.2f}")


def plot(rows: list[Row], path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is not installed; rerun without --plot or install matplotlib") from exc

    path.parent.mkdir(parents=True, exist_ok=True)
    n_values = sorted({r.n for r in rows})
    fig, ax = plt.subplots(figsize=(11, 6))
    markers = {"prime": "o", "prime_power": "s", "composite": "."}

    for n in n_values:
        n_rows = [r for r in rows if r.n == n]
        for k in sorted({r.kind for r in n_rows}):
            items = [r for r in n_rows if r.kind == k]
            ax.scatter(
                [r.q for r in items],
                [r.ratio for r in items],
                label=f"{k}, n={n}",
                marker=markers.get(k, "."),
                s=28 if k != "composite" else 14,
                alpha=0.75,
            )

    ax.set_title("Carmichael lambda ratio scan: lambda(q) / log(n)")
    ax.set_xlabel("q")
    ax.set_ylabel("lambda(q) / log(n)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize="small", ncols=2)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    print(f"\nWrote plot: {path}")


def parse_n_values(raw: str) -> list[int]:
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Explore R(q,n)=lambda(q)/log(n).")
    parser.add_argument("--q-min", type=int, default=2)
    parser.add_argument("--q-max", type=int, default=50)
    parser.add_argument("--n", type=parse_n_values, default=[1000], help="comma-separated n values, e.g. 100,1000,10000")
    parser.add_argument("--method", choices=["fast", "brute"], default="fast")
    parser.add_argument("--csv", type=Path, help="optional CSV output path")
    parser.add_argument("--plot", type=Path, help="optional PNG plot output path; requires matplotlib")
    parser.add_argument("--top", type=int, default=10, help="number of ratio spikes to summarize")
    parser.add_argument("--limit", type=int, help="limit printed table rows")
    args = parser.parse_args()

    rows = scan(args.q_min, args.q_max, args.n, method=args.method)
    print_table(rows, limit=args.limit)
    summarize(rows, top=args.top)

    if args.csv:
        write_csv(rows, args.csv)
        print(f"\nWrote CSV: {args.csv}")
    if args.plot:
        plot(rows, args.plot)


if __name__ == "__main__":
    main()
