---
title: "Multiplicative Collapse in $(\\mathbb{Z}/n\\mathbb{Z})^*$: Wedge Envelopes and Synchronization Growth"
author: "lambda_ratio_explorer"
geometry: margin=1in
---

# Multiplicative Collapse in $(\mathbb{Z}/n\mathbb{Z})^*$

> Factorization splits the multiplicative group into independent cyclic
> components. $\varphi(n)$ counts the total elements; $\lambda(n)$ is the
> maximum order of any element. Their ratio measures how badly those
> components fail to synchronize into a single cycle, and that obstruction
> is governed by the gcd of the component cycle lengths.

This note states two elementary identities that organize the visualizations
in this repository: a wedge envelope for $\lambda$, and a multiplicative
propagation law for the collapse index $C = \varphi/\lambda$. Both follow
from standard properties of $\varphi$ and $\lambda$ in two or three lines.
The contribution is the framing, not new mathematics. The honesty note at
the end says exactly that.

## 1. Setup

For $n \ge 1$ let $\varphi(n) = \lvert (\mathbb{Z}/n\mathbb{Z})^* \rvert$ be
Euler's totient and let $\lambda(n)$ denote the Carmichael function: the
exponent of the multiplicative group, that is, the smallest positive integer
$m$ with $a^m \equiv 1 \pmod{n}$ for every $a$ coprime to $n$. Define the
**collapse index**

$$
C(n) \;=\; \frac{\varphi(n)}{\lambda(n)}.
$$

Then $C(n)$ is a positive integer (since $\lambda(n) \mid \varphi(n)$) and
$C(n) = 1$ if and only if $(\mathbb{Z}/n\mathbb{Z})^*$ is cyclic. For
semiprimes $n = pq$ with distinct odd primes $p, q$,

$$
C(pq) \;=\; \gcd(p - 1,\; q - 1),
$$

which is exactly the cryptographic-strength dial in RSA-style moduli.

## 2. Standard lemmas

**Lemma 1 (multiplicativity of $\varphi$).** If $\gcd(a, b) = 1$, then
$\varphi(ab) = \varphi(a)\varphi(b)$.

**Lemma 2 (Carmichael under coprime products).** If $\gcd(a, b) = 1$, then
$\lambda(ab) = \mathrm{lcm}(\lambda(a), \lambda(b))$.

Both are standard. Lemma 2 follows from the Chinese Remainder Theorem
isomorphism $(\mathbb{Z}/ab\mathbb{Z})^* \cong (\mathbb{Z}/a\mathbb{Z})^*
\times (\mathbb{Z}/b\mathbb{Z})^*$ and the fact that the exponent of a
direct product of finite abelian groups is the lcm of their exponents.

We will also use the elementary identity

$$
\mathrm{lcm}(x, y) \;=\; \frac{x \cdot y}{\gcd(x, y)}.
$$

## 3. Wedge envelope

**Theorem 1 (Wedge envelope).** Let $q = k \cdot p$ with $p$ prime and
$p \nmid k$. Then

$$
\lambda(q) \;=\; \mathrm{lcm}\bigl(\lambda(k),\; p - 1\bigr).
$$

Two consequences are useful:

1. **Upper bound.** $\lambda(q) \le \lambda(k) \cdot (p - 1)$, with equality iff $\gcd(\lambda(k), p - 1) = 1$.
2. **Envelope line.** If $\lambda(k) \mid (p - 1)$, then $\lambda(q) = p - 1$. When additionally $p \ge k$, this gives the linear envelope $\lambda(q) = p - 1 \le q/k - 1$.

**Proof.** $\lambda(q) = \mathrm{lcm}(\lambda(k), p - 1)$ is Lemma 2.
For (1), the standard inequality $\mathrm{lcm}(x, y) \le xy$ becomes equality iff
$\gcd(x, y) = 1$. For (2), if $\lambda(k) \mid (p - 1)$ then $\mathrm{lcm}(\lambda(k), p - 1) = p - 1$;
the bound $p - 1 \le q/k - 1$ rearranges to $kp \ge k(p - 1) + k$, i.e. trivially $q \ge q$. $\blacksquare$

**Where the envelope is tight.** The case $k = 2$ is automatic: $\lambda(2) = 1$
divides every $p - 1$, so for every odd prime $p$ the semiprime $q = 2p$ sits
exactly on the line $\lambda = q/2 - 1$. The case $k = 3$ is also automatic
for $p > 2$, since $\lambda(3) = 2$ divides $p - 1$ for any odd prime. For
$k = 5$ ($\lambda(5) = 4$), the envelope is hit only by primes $p \equiv 1 \pmod 4$,
which is half of them; the others land below.

These are the diagonal rays in `wedge_envelopes.png`. Each line
$\lambda = q/k - 1$ is the upper boundary of the points whose smallest prime
factor is $k$, populated exactly by those semiprimes for which the
divisibility $\lambda(k) \mid (p - 1)$ is satisfied.

## 4. Collapse propagation

**Theorem 2 (Collapse propagation).** Let $q = k \cdot p$ with $p$ prime
and $p \nmid k$. Then

$$
\boxed{\;C(q) \;=\; C(k) \cdot \gcd\!\bigl(\lambda(k),\; p - 1\bigr).\;}
$$

**Proof.** By Lemma 1, $\varphi(q) = \varphi(k) \cdot (p - 1)$, since
$\gcd(k, p) = 1$. By Lemma 2, $\lambda(q) = \mathrm{lcm}(\lambda(k), p - 1)$.
Then

$$
C(q) \;=\; \frac{\varphi(q)}{\lambda(q)}
\;=\; \frac{\varphi(k)(p - 1)}{\mathrm{lcm}(\lambda(k), p - 1)}
\;=\; \frac{\varphi(k)(p - 1) \cdot \gcd(\lambda(k), p - 1)}
{\lambda(k)(p - 1)}
\;=\; \frac{\varphi(k)}{\lambda(k)} \cdot \gcd(\lambda(k), p - 1)
\;=\; C(k) \cdot \gcd(\lambda(k), p - 1). \qquad \blacksquare
$$

**Interpretation.** Adding a new prime factor $p$ multiplies the collapse
index by exactly the amount its cycle length $p - 1$ overlaps with the
existing exponent $\lambda(k)$. There is no smooth shrinkage: collapse
moves in clean integer multiplicative jumps, and each jump is a gcd.

**Iterating.** Starting from $k_0 = 1$ with $C(1) = \lambda(1) = 1$, applying
Theorem 2 to a sequence of distinct primes $p_1, p_2, \ldots, p_r$ gives

$$
C(p_1 p_2 \cdots p_r) \;=\; \prod_{i=1}^{r} \gcd\!\bigl(\lambda(p_1 \cdots p_{i-1}),\; p_i - 1\bigr),
$$

with the convention $\lambda(1) = 1$ at the first step. The order of the
primes does not change the final $C$ since both $\varphi$ and $\lambda$
are determined by the multiset of factors.

## 5. Worked examples

### 5.1 $1365 = 3 \cdot 5 \cdot 7 \cdot 13$

| step | $k$  | $p$ | $\lambda(k)$ | $\gcd(\lambda(k), p-1)$ | $kp$  | $C(kp)$ |
| ---- | ---- | --- | ------------ | ----------------------- | ----- | ------- |
| 1    | 1    | 3   | 1            | 1                       | 3     | 1       |
| 2    | 3    | 5   | 2            | 2                       | 15    | 2       |
| 3    | 15   | 7   | 4            | 2                       | 105   | 4       |
| 4    | 105  | 13  | 12           | **12**                  | 1365  | **48**  |

The dramatic last step is because $\lambda(105) = 12$ already contains a
full $4 \cdot 3$ structure, and $13 - 1 = 12$ matches it exactly.

### 5.2 $1729 = 7 \cdot 13 \cdot 19$ (Carmichael)

| step | $k$  | $p$ | $\lambda(k)$ | $\gcd(\lambda(k), p-1)$ | $kp$  | $C(kp)$ |
| ---- | ---- | --- | ------------ | ----------------------- | ----- | ------- |
| 1    | 1    | 7   | 1            | 1                       | 7     | 1       |
| 2    | 7    | 13  | 6            | 6                       | 91    | 6       |
| 3    | 91   | 19  | 12           | 6                       | 1729  | **36**  |

That $1729$ is a Carmichael number is a separate condition
$\lambda(n) \mid n - 1$, but the high collapse index is a *consequence* of
the same shared structure: $7 - 1 = 6$, $13 - 1 = 12$, $19 - 1 = 18$ all
share the factors of $6$.

### 5.3 Two semiprimes of similar magnitude

| $n$  | factors        | $\gcd(p - 1, q - 1)$ | $C(n)$ |
| ---- | -------------- | --------------------- | ------ |
| 253  | $11 \cdot 23$  | $\gcd(10, 22) = 2$    | 2      |
| 91   | $7 \cdot 13$   | $\gcd(6, 12) = 6$     | 6      |

The two have nearly the same magnitude but very different group structure.
The clean shape is $253$: the two primes' totients share only the inevitable
factor of $2$. The overlapping shape is $91$: the totients share both 2 and 3,
giving a threefold-larger collapse.

## 6. Distribution of collapse for semiprimes

Theorems 1 and 2 are exact identities. The next question is *how big*
$C(pq)$ tends to be when $p, q$ are random primes. The answer is governed
by Dirichlet's theorem on primes in arithmetic progressions.

**Theorem C (Dirichlet density of $\gcd(p-1, q-1)$).** Let $p, q$ be
distinct odd primes drawn uniformly from primes $\le X$, and let $\ell$
be a prime. Then

$$
\Pr\!\bigl(\ell \mid \gcd(p - 1,\; q - 1)\bigr) \;\xrightarrow[X \to \infty]{}\; \frac{1}{(\ell - 1)^2}.
$$

The case $\ell = 2$ is trivial: every odd prime $p$ has $p - 1$ even, so
the probability is $1$ (and the formula gives $1/(2-1)^2 = 1$ as a
boundary value).

**Justification.** By Dirichlet's theorem, for any prime $\ell$ the primes
$p$ with $p \equiv 1 \pmod{\ell}$ have natural density $1/(\ell - 1)$
among all primes. The two prime draws are asymptotically independent, so
the joint event $\ell \mid (p - 1)$ and $\ell \mid (q - 1)$ has density
$1/(\ell - 1)^2$. The same argument with prime powers $\ell^k$ gives
$\Pr(\ell^k \mid \gcd) \to 1/\bigl(\ell^{k-1}(\ell - 1)\bigr)^2$ (for $\ell$ odd).

**Corollary (asymptotic expectation).** Using the divisor identity
$\gcd(a, b) = \sum_{d \ge 1} \varphi(d) \, [d \mid a]\,[d \mid b]$ together
with the Dirichlet rate $\Pr(d \mid p - 1) \to 1/\varphi(d)$,

$$
\mathbb{E}\bigl[\gcd(p - 1,\; q - 1)\bigr] \;\approx\; \sum_{d \le X} \frac{1}{\varphi(d)} \;\sim\; A \, \log X,
$$

where $A = \dfrac{315\,\zeta(3)}{2\pi^4} \approx 1.9436$ is the standard
mean-totient-reciprocal constant. The expectation is unbounded but grows
only logarithmically with the prime cutoff. The dominant terms come from
small primes: $d = 1, 2, 3, 4, 6$ alone contribute well over half of the
sum at any $X$, which is why the *shape* of Theorem C, not just its
asymptotic, is dominated by $\ell \in \{2, 3\}$.

**Mechanism of the divergence.** The growth of the sum is not driven by
any single dominant divisor: $\sum_{d \le X} 1/\varphi(d) \sim A \log X$
accumulates slowly because integers with small totients (those with many
small prime factors) become rare, but they remain frequent enough at
every scale to keep the partial sum unbounded. The divergence arises
from the slow density of integers with small totients, not from any
single dominant divisor.

**Why $1365 = 3 \cdot 5 \cdot 7 \cdot 13$ is at the top.** Theorem C
predicts that the most-collapsed semiprimes are products of *small*
primes whose $(p_i - 1)$ values share heavy small-prime factors. For
$1365$, the four primes $3, 5, 7, 13$ contribute totients $2, 4, 6, 12$,
all rich in factors of $2$ and $3$. Each step of the propagation trace
in §5.1 picks up the maximum possible $\gcd$ at that step. Theorem C is
the statistical version of the same statement: if you sample primes
uniformly, the expected gcd is dominated by exactly the kinds of small
shared structure that $1365$ exhibits in concentrated form.

The empirical verification of Theorem C, including a divisibility-rate
table for $\ell \in \{2, 3, 5, 7, 11, 13, 17, 19, 23\}$, is in
`gcd_distribution_theory.py`. The script prints empirical vs predicted
rates side by side and writes `gcd_distribution.png`.

## 7. Cryptographic interpretation

For RSA, the modulus $n = pq$ is chosen with $p, q$ large distinct primes.
The decryption exponent $d$ must satisfy $ed \equiv 1 \pmod{\lambda(n)}$
rather than $\pmod{\varphi(n)}$, and $\lambda(n) = \varphi(n) / \gcd(p-1, q-1)$.
A large $\gcd(p-1, q-1)$ shrinks $\lambda$ relative to $\varphi$, shortens
all element cycles, and accelerates Pohlig-Hellman style attacks on the
multiplicative group.

In the language of Theorem 2:

> RSA key generation tries to choose the second prime so that
> $\gcd(\lambda(k), p - 1)$ is as small as possible at the moment of joining.

The "safe prime" condition that $(p - 1)/2$ is also prime is one way to
guarantee this: it forces $\gcd(p - 1, q - 1) \in \{1, 2\}$ for any other
safe prime $q$, since $p - 1$ has only the factors $1, 2, (p-1)/2, p-1$.

The heatmap panel of `group_structure.png` is exactly this map. Bright
cells are pairs $(p, q)$ with high $\gcd(p - 1, q - 1)$; dark cells are
the cryptographically clean pairs.

## 8. Why this matters

The three theorems are not independent observations. They are the same
multiplicative law shown at three scales:

- **RSA parameter hygiene.** Safe-prime selection is the operational
  version of "minimize $\gcd(\lambda(k), p - 1)$ at the moment of joining."
  Theorem 2 is the algebraic statement of that goal; Theorem C explains
  why the goal is hard: although most pairs of primes share only small
  factors ($2$, $3$, $5$), the cumulative effect of all possible shared
  divisors causes the expected value of $\gcd(p - 1, q - 1)$ to grow
  logarithmically with prime size. Safe-prime selection prevents this
  accumulation by structurally restricting the divisor lattice, forcing
  the gcd into $\{1, 2\}$ regardless of $X$.
- **Carmichael behavior.** Carmichael numbers are the extreme case of
  the same shared-structure phenomenon. The condition
  $\lambda(n) \mid n - 1$ requires precisely the kind of heavy
  $\gcd$-overlap that maximizes $C(n)$. Korselt's criterion is a
  consequence; the gcd story is the cause.
- **Visual intuition.** The wedge envelopes (`wedge_envelopes.png`),
  the collapse bands (`group_structure.png`), the element-order
  histograms (`order_distributions.png`), and the propagation steps
  (`propagation.png`, `gcd_distribution.png`) are all the same law
  rendered at different scales. Once you see one, you see them all.

## 9. Connection to the figures in this repository

| File                       | What it shows                                                                  |
| -------------------------- | ------------------------------------------------------------------------------ |
| `wedge_envelopes.png`      | Theorem 1 made visible: lines $\lambda = q/\ell - 1$ for each smallest factor $\ell$. |
| `group_structure.png`      | Distribution of $C(q)$, the gcd-heatmap, and Carmichael overlays.              |
| `order_distributions.png`  | The parallel-cycle picture inside specific $n$, with explicit invariant factors. |
| `propagation.png`          | Theorem 2 made visible: stacked $\log_2 \gcd$ contributions per prime, plus the empirical density of $\gcd(p - 1, q - 1)$. |
| `gcd_distribution.png`     | Theorem C made visible: empirical vs Dirichlet-predicted divisibility rates, plus the $1365$ tie-back. |

## 10. Honest scope

Theorems 1 and 2 are elementary three-line consequences of the standard
multiplicativity of $\varphi$ and $\lambda$. They are well known in
substance to anyone fluent with these functions, and they would not pass
peer review at a research journal as new mathematics. Theorem C is a
direct application of Dirichlet's theorem on primes in arithmetic
progressions and is similarly well known in substance; the contribution
of this note is to state it in the same vocabulary as Theorems 1 and 2
and to verify it empirically against the same examples.

What this note offers is a clean *framing*: named theorems that organize
the wedge envelopes, the collapse bands, the synchronization model, and
the element-order distributions into one multiplicative law. The repository
treats these results as teaching artifacts and a structural spine, not as
research. Scope is exactly what is stated and proved here.

The value, in one sentence: collapse propagates multiplicatively, and the
multiplier is a gcd you can see in the picture.
