# lambda_ratio_explorer

A small numerical and visual explorer for the multiplicative group `(Z/nZ)*`.
The central object is the **collapse index**

```text
C(n) = phi(n) / lambda(n)
```

where `phi` is Euler's totient and `lambda` is the Carmichael function.

`C(n)` counts how many cyclic factors the multiplicative group decomposes
into beyond its largest one. It is exactly **1** when the group is cyclic
(every prime, every odd prime power, plus a handful of small special cases)
and grows with structural collapse.

For semiprimes `n = p*q` it has a clean closed form:

```text
C(p*q) = gcd(p - 1, q - 1)
```

This is precisely the quantity that controls cycle overlap in `(Z/nZ)*` and
is the cryptographic-strength dial in RSA-style moduli.

## Install

```bash
python -m venv .venv
.venv\Scripts\activate          # PowerShell on Windows
# source .venv/bin/activate     # bash / zsh
pip install -r requirements.txt
```

## Run

### Core scanner

`lambda_ratio_explorer.py` is the CLI. It computes `phi`, `lambda`, and
`C` for each `q` in a range and writes a table, optional CSV, and
optional scatter plot.

```bash
python lambda_ratio_explorer.py --q-max 200 --n 1000 \
    --csv runs/scan.csv --plot runs/scan.png
```

### Group-structure analysis (`group_structure.py`)

The teaching centerpiece. Generates a 4-panel figure showing collapse
behavior across q, the prime-pair heatmap of `gcd(p-1, q-1)`, the
distribution of `C` by structural kind, and a side-by-side comparison
of two collapse measures.

```bash
python group_structure.py
# -> group_structure.png
```

### Wedge envelopes (`wedge_envelopes.py`)

A focused plot showing that the upper edges of the `lambda(q)` wedge
are **exact** algebraic lines `q/k - 1`, indexed by the smallest prime
factor `k` of `q`. Useful as an algebraic warm-up before the collapse
discussion.

```bash
python wedge_envelopes.py
# -> wedge_envelopes.png
```

## What the data shows

Running `group_structure.py` over `q in [2, 2000]`:

| kind         | count | mean C | median C | max C | fraction with C=1 |
| ------------ | ----- | ------ | -------- | ----- | ----------------- |
| prime        | 303   | 1.00   | 1        | 1     | 1.000             |
| prime_power  | 30    | 1.27   | 1        | 2     | 0.733             |
| composite    | 1666  | 4.96   | 4        | 48    | 0.110             |

Primes are exactly cyclic. Composites collapse, and the collapse is
quantized in integer tiers. The most-collapsed values in the range are
products of small primes whose `(p_i - 1)` shares many common factors
(e.g. `1365 = 3 * 5 * 7 * 13` with `C = 48`).

The Hardy-Ramanujan number `1729 = 7 * 13 * 19` shows up near the top:
it is also a Carmichael number, and large `C` is part of why.

## Cryptographic interpretation

For `n = p*q` (an RSA-shaped modulus), the order of an arbitrary unit
divides `lambda(n) = lcm(p-1, q-1)`, not `phi(n) = (p-1)(q-1)`. The
gap between them is `C(n) = gcd(p-1, q-1)`. Choosing `p`, `q` so that
`C(n)` is small (ideally 2) is part of what makes a modulus
cryptographically clean.

The heatmap panel of `group_structure.png` is exactly this map: bright
cells are prime pairs to avoid, dark cells are pairs whose totients
share little.

## Files

| File                       | Purpose                                           |
| -------------------------- | ------------------------------------------------- |
| `lambda_ratio_explorer.py` | Core library + CLI scanner                        |
| `group_structure.py`       | Centerpiece 4-panel analysis of `C(n)`            |
| `wedge_envelopes.py`       | Algebraic envelope visualization                  |
| `requirements.txt`         | `matplotlib` (pulls in numpy)                     |
| `lambda_ratio_scan_*.csv`  | Sample CSV outputs from earlier scanner runs      |

## Notes

- Factorization is trial division. Practical up to `q` of order 10^5.
  Replace with Pollard rho if you want to push further.
- The legacy ratio `lambda(q) / log(n)` is retained in `Row` and the CLI
  for backward compatibility, but `log(n)` is just a constant scalar
  and does not enter the structural story. The interesting metrics are
  `lambda(q)`, `phi(q)`, and `C(q)`.
