"""Tests for the core library in lambda_ratio_explorer.py.

The brute-force Carmichael implementation serves as the oracle for the
fast one, and the THEOREM.md worked examples serve as fixtures for the
collapse propagation theorem.
"""

import sys
from collections import Counter
from math import gcd
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from lambda_ratio_explorer import (
    carmichael_lambda,
    collapse_index,
    collapse_propagation_trace,
    collapse_step,
    divisors,
    element_orders,
    euler_totient,
    factorize,
    fracture_count,
    invariant_factors,
    is_carmichael,
    is_prime,
    kind,
)


PRIMES_TO_50 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# Carmichael numbers below 10^4 (OEIS A002997).
CARMICHAEL_BELOW_10K = [561, 1105, 1729, 2465, 2821, 6601, 8911]


class TestFactorize:
    @pytest.mark.parametrize("n", range(1, 500))
    def test_reconstructs_n(self, n):
        product = 1
        for p, e in factorize(n).items():
            assert is_prime(p)
            product *= p ** e
        assert product == n


class TestIsPrime:
    def test_matches_sieve_below_10000(self):
        limit = 10000
        sieve = [True] * (limit + 1)
        sieve[0] = sieve[1] = False
        for i in range(2, int(limit ** 0.5) + 1):
            if sieve[i]:
                for j in range(i * i, limit + 1, i):
                    sieve[j] = False
        for n in range(limit + 1):
            assert is_prime(n) == sieve[n]

    def test_large_primes_and_composites(self):
        assert is_prime(2 ** 61 - 1)  # Mersenne prime
        assert is_prime(1000000007)
        assert not is_prime(1000000007 * 1000000009)
        # strong pseudoprime to base 2, composite
        assert not is_prime(3215031751)
        # psi_12: smallest strong pseudoprime to all 12 prime bases <= 37;
        # base 41 in the witness set must catch it
        assert not is_prime(318665857834031151167461)


class TestFactorizeLarge:
    def test_fermat_number_f6(self):
        # F6 = 2^64 + 1 = 274177 * 67280421310721, both prime
        assert factorize(2 ** 64 + 1) == {274177: 1, 67280421310721: 1}

    def test_18_digit_semiprime(self):
        p, q = 1000000007, 1000000009
        assert factorize(p * q) == {p: 1, q: 1}

    def test_large_prime_power(self):
        p = 1000003
        assert factorize(p ** 3) == {p: 3}

    def test_semiprime_collapse_identity_large(self):
        p, q = 1000000007, 1000000009
        assert collapse_index(p * q) == gcd(p - 1, q - 1)


class TestCarmichaelLambda:
    @pytest.mark.parametrize("n", range(1, 151))
    def test_fast_matches_bruteforce(self, n):
        assert carmichael_lambda(n) == carmichael_lambda(n, method="brute")

    def test_powers_of_two(self):
        # lambda(2) = 1, lambda(4) = 2, lambda(2^k) = 2^(k-2) for k >= 3
        assert carmichael_lambda(2) == 1
        assert carmichael_lambda(4) == 2
        assert carmichael_lambda(8) == 2
        assert carmichael_lambda(16) == 4
        assert carmichael_lambda(1024) == 256

    @pytest.mark.parametrize("n", range(2, 300))
    def test_divides_phi(self, n):
        assert euler_totient(n) % carmichael_lambda(n) == 0


class TestCollapseIndex:
    def test_semiprime_gcd_identity(self):
        odd_primes = [p for p in PRIMES_TO_50 if p > 2]
        for i, p in enumerate(odd_primes):
            for q in odd_primes[i + 1:]:
                assert collapse_index(p * q) == gcd(p - 1, q - 1)

    @pytest.mark.parametrize("p", PRIMES_TO_50)
    def test_primes_are_cyclic(self, p):
        assert collapse_index(p) == 1

    def test_known_extremes(self):
        assert collapse_index(1365) == 48  # 3 * 5 * 7 * 13
        assert collapse_index(1729) == 36  # 7 * 13 * 19


class TestInvariantFactors:
    @pytest.mark.parametrize("n", range(2, 300))
    def test_structure(self, n):
        factors = invariant_factors(n)
        assert factors == sorted(factors)
        # divisibility chain d_1 | d_2 | ... | d_k
        for a, b in zip(factors, factors[1:]):
            assert b % a == 0
        product = 1
        for d in factors:
            product *= d
        assert product == euler_totient(n)
        assert factors[-1] == carmichael_lambda(n)
        assert fracture_count(n) == len(factors)

    @pytest.mark.parametrize("n", [7, 8, 15, 16, 24, 35, 63, 91, 105])
    def test_matches_element_orders(self, n):
        """The order histogram of the abstract product of cyclic groups
        must match the actual element orders in (Z/nZ)*.

        In a cyclic group of order d there are phi(e) elements of order e
        for each e | d; in a product, orders combine by lcm.
        """
        expected = Counter({1: 1})
        for d in invariant_factors(n):
            combined = Counter()
            for order, count in expected.items():
                for e in divisors(d):
                    combined[lcm_pair(order, e)] += count * euler_totient(e)
            expected = combined
        actual = Counter(element_orders(n).values())
        assert actual == expected


def lcm_pair(a, b):
    return a * b // gcd(a, b)


class TestCollapsePropagation:
    def test_theorem_worked_example_1365(self):
        trace = collapse_propagation_trace([3, 5, 7, 13])
        assert [s["C_kp"] for s in trace] == [1, 2, 4, 48]
        assert trace[-1]["kp"] == 1365

    def test_theorem_worked_example_1729(self):
        trace = collapse_propagation_trace([7, 13, 19])
        assert [s["C_kp"] for s in trace] == [1, 6, 36]
        assert trace[-1]["kp"] == 1729

    def test_order_independence(self):
        for perm in ([13, 7, 5, 3], [5, 13, 3, 7]):
            assert collapse_propagation_trace(perm)[-1]["C_kp"] == 48

    def test_rejects_duplicate_prime(self):
        with pytest.raises(ValueError):
            collapse_propagation_trace([3, 3])

    def test_rejects_composite(self):
        with pytest.raises(ValueError):
            collapse_step(1, 6)

    def test_rejects_shared_factor(self):
        with pytest.raises(ValueError):
            collapse_step(15, 5)


class TestIsCarmichael:
    def test_known_carmichael_numbers(self):
        found = [n for n in range(2, 10000) if is_carmichael(n)]
        assert found == CARMICHAEL_BELOW_10K

    @pytest.mark.parametrize("n", range(2, 3000))
    def test_matches_lambda_criterion(self, n):
        """Carmichael numbers are exactly the composites with lambda(n) | n - 1
        (and more than one prime factor, i.e. squarefree composites)."""
        expected = (
            not is_prime(n)
            and len(factorize(n)) > 1
            and all(e == 1 for e in factorize(n).values())
            and (n - 1) % carmichael_lambda(n) == 0
        )
        assert is_carmichael(n) == expected


class TestCollapseAtScaleScan:
    def test_sieve_scan_matches_library(self):
        pytest.importorskip("matplotlib")
        from collapse_at_scale import scan

        C, kind_code = scan(2000)
        kind_names = {1: "prime", 2: "prime_power", 3: "composite"}
        for n in range(2, 2001):
            assert C[n] == collapse_index(n)
            assert kind_names[kind_code[n]] == kind(n)


class TestKind:
    def test_classification(self):
        assert kind(7) == "prime"
        assert kind(8) == "prime_power"
        assert kind(9) == "prime_power"
        assert kind(12) == "composite"
