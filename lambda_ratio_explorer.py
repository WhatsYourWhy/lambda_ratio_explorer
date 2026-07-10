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
import random
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


# Factors below this bound are stripped by trial division; anything left
# is handled by Miller-Rabin + Pollard rho.
_TRIAL_DIVISION_LIMIT = 10_000


def _pollard_rho(n: int) -> int:
    """Return a nontrivial factor of odd composite n (Floyd cycle detection)."""
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d


def _factor_hard(n: int) -> list[int]:
    """Prime factors (with multiplicity) of n, which has no factor below
    _TRIAL_DIVISION_LIMIT."""
    if n == 1:
        return []
    if is_prime(n):
        return [n]
    d = _pollard_rho(n)
    return _factor_hard(d) + _factor_hard(n // d)


def factorize(n: int) -> dict[int, int]:
    """Prime factorization: trial division for small factors, then
    Pollard rho for the remaining cofactor. Practical well beyond the
    old pure-trial-division limit (e.g. 18-digit semiprimes)."""
    if n < 1:
        raise ValueError("n must be positive")

    factors: dict[int, int] = {}
    d = 2
    while d * d <= n and d <= _TRIAL_DIVISION_LIMIT:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2  # 2, then odd candidates only
    for p in _factor_hard(n):
        factors[p] = factors.get(p, 0) + 1
    return dict(sorted(factors.items()))


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


def divisors(n: int) -> list[int]:
    """Sorted list of all positive divisors of n."""
    if n < 1:
        raise ValueError("n must be positive")
    if n == 1:
        return [1]
    result = [1]
    for p, e in factorize(n).items():
        new_divs: list[int] = []
        power = 1
        for _ in range(e + 1):
            for d in result:
                new_divs.append(d * power)
            power *= p
        result = new_divs
    return sorted(result)


def _cyclic_factors_of_unit_group_at_prime_power(p: int, k: int) -> list[int]:
    """Cyclic factor orders of (Z/p^k Z)*.

    For odd p, (Z/p^k Z)* is cyclic of order phi(p^k).
    For p = 2: trivial at k=1, cyclic of order 2 at k=2,
               and Z/2 x Z/2^(k-2) for k >= 3.
    """
    if p == 2:
        if k == 1:
            return []
        if k == 2:
            return [2]
        return [2, 2 ** (k - 2)]
    return [(p - 1) * p ** (k - 1)]


def invariant_factors(n: int) -> list[int]:
    """Invariant factor decomposition of (Z/nZ)*.

    Returns d_1 | d_2 | ... | d_k with d_1 <= ... <= d_k = lambda(n) and
    d_1 * d_2 * ... * d_k = phi(n). The length k is the fracture count:
    the number of cyclic components in the decomposition.
    """
    if n < 1:
        raise ValueError("n must be positive")
    if n <= 2:
        return [1]

    cyclic_orders: list[int] = []
    for p, k in factorize(n).items():
        cyclic_orders.extend(_cyclic_factors_of_unit_group_at_prime_power(p, k))
    if not cyclic_orders:
        return [1]

    prime_to_vals: dict[int, list[int]] = {}
    for order in cyclic_orders:
        for q, e in factorize(order).items():
            prime_to_vals.setdefault(q, []).append(e)

    if not prime_to_vals:
        return [1]

    k_max = max(len(v) for v in prime_to_vals.values())
    for q in prime_to_vals:
        prime_to_vals[q].sort(reverse=True)
        while len(prime_to_vals[q]) < k_max:
            prime_to_vals[q].append(0)

    factors: list[int] = []
    for j in range(k_max):
        d = 1
        for q, vals in prime_to_vals.items():
            d *= q ** vals[j]
        if d > 1:
            factors.append(d)

    factors.sort()
    return factors if factors else [1]


def fracture_count(n: int) -> int:
    """Number of cyclic components in the invariant factor decomposition."""
    return len(invariant_factors(n))


def element_orders(n: int) -> dict[int, int]:
    """Map each unit a in (Z/nZ)* to its multiplicative order.

    Computed by testing divisors of lambda(n) in ascending order, since
    every element order divides lambda(n).
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    lam = carmichael_lambda(n)
    divs = divisors(lam)
    orders: dict[int, int] = {}
    for a in range(1, n):
        if gcd(a, n) != 1:
            continue
        for d in divs:
            if pow(a, d, n) == 1:
                orders[a] = d
                break
    return orders


def collapse_step(k: int, p: int) -> dict:
    """One step of the collapse propagation theorem.

    For prime p with p not dividing k, the multiplicative collapse index
    satisfies the identity

        C(k * p) = C(k) * gcd(lambda(k), p - 1)

    Returns a dict with all quantities relevant to the step:
        k, p, phi_k, lambda_k, C_k,
        gcd_term, kp, phi_kp, lambda_kp, C_kp.

    Raises ValueError if p is not prime, if p divides k, or if k < 1.
    The identity is asserted at runtime so this function doubles as a test.
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    if not is_prime(p):
        raise ValueError(f"p={p} is not prime")
    if k % p == 0:
        raise ValueError(f"p={p} divides k={k}; theorem requires gcd(k, p) = 1")

    phi_k = euler_totient(k)
    lambda_k = carmichael_lambda(k)
    C_k = phi_k // lambda_k

    gcd_term = gcd(lambda_k, p - 1)
    kp = k * p
    phi_kp = euler_totient(kp)
    lambda_kp = carmichael_lambda(kp)
    C_kp = phi_kp // lambda_kp

    predicted = C_k * gcd_term
    assert C_kp == predicted, (
        f"propagation identity violated: C({kp}) = {C_kp} but "
        f"C({k}) * gcd(lambda({k}), {p}-1) = {C_k} * {gcd_term} = {predicted}"
    )

    return {
        "k": k,
        "p": p,
        "phi_k": phi_k,
        "lambda_k": lambda_k,
        "C_k": C_k,
        "gcd_term": gcd_term,
        "kp": kp,
        "phi_kp": phi_kp,
        "lambda_kp": lambda_kp,
        "C_kp": C_kp,
    }


def collapse_propagation_trace(primes: list[int]) -> list[dict]:
    """Iterate collapse_step starting at k = 1 across the given primes in order.

    Returns the full trace as a list of step dicts. Each successive step uses
    the previous step's k*p as the new k. The primes must be distinct and the
    order matters only in that earlier values appear in lambda_k of later
    steps; the final C is independent of order.
    """
    if not primes:
        return []
    seen: set[int] = set()
    trace: list[dict] = []
    k = 1
    for p in primes:
        if p in seen:
            raise ValueError(f"duplicate prime {p} in propagation trace")
        seen.add(p)
        step = collapse_step(k, p)
        trace.append(step)
        k = step["kp"]
    return trace


def is_carmichael(n: int) -> bool:
    """Korselt's criterion: n is a Carmichael number iff n is composite,
    squarefree, and (p - 1) divides (n - 1) for every prime p dividing n.

    Carmichael numbers are exactly the composites with lambda(n) | (n - 1),
    so they pass the Fermat primality test for every coprime base. The
    smallest is 561 = 3 * 11 * 17.
    """
    if n < 561 or is_prime(n):
        return False
    factors = factorize(n)
    if len(factors) < 3:
        return False
    if any(e > 1 for e in factors.values()):
        return False
    return all((n - 1) % (p - 1) == 0 for p in factors)


# Deterministic Miller-Rabin witness set: correct for all n < 3.3 * 10^24.
_MILLER_RABIN_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin, exact for all n < 3.3 * 10^24."""
    if n < 2:
        return False
    for p in _MILLER_RABIN_WITNESSES:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in _MILLER_RABIN_WITNESSES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
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
