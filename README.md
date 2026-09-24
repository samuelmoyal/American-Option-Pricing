# American-Option-Pricing


For the theory, please refer to the [project report](Report.pdf) (March 2024). The numerical results of the
report were produced by a buggy first version of the code; they are superseded by the notebooks and by the
[Results](#results) section below (see [Corrections](#corrections)).

# Summary


This repository contains a comprehensive analysis and implementation of various methods for pricing American options. The project is designed to provide a clear understanding of American option pricing mechanisms and to implement these using mathematical and computational techniques.

## Table of Contents

- [Introduction](#introduction)
- [Methods](#methods)
- [Usage](#usage)
- [Results](#results)
- [Corrections](#corrections)



## Introduction

American options are financial derivatives that give the holder the right to exercise the option at any time before or at the expiration date. This flexibility makes American options more complex to price than European options, as they require consideration of optimal exercise strategies.

This project explores different methods to price American options, comparing their efficiency, accuracy, and computational demands.

## Methods

The repository implements and analyzes the following pricing methods:

1. **Binomial Trees** (Cox-Ross-Rubinstein): European, American and Bermudan versions, used as the benchmark.
2. **Quantization Tree** (Bally-Pagès): each date is split into equally likely cells, with exact transition probabilities between cells.
3. **Least Squares Monte Carlo** (Longstaff-Schwartz): regression of the continuation value on Laguerre polynomials of S/K, on in-the-money paths.

All pricers live in [`pricing.py`](pricing.py); the notebooks run the studies.

Each method includes explanations, code implementations, and analysis of its advantages and disadvantages.

## Usage

Requirements: `numpy`, `scipy`, `matplotlib` (and Jupyter to run the notebooks).

- **Binomial Trees**: Run the notebook `Binomial_tree.ipynb`
- **Discretization Method**: Run the notebook `DiscretisationVF.ipynb`
- **Monte Carlo Simulation**: Run the notebook `Monte_carlo.ipynb`

Each notebook includes configurable parameters for option characteristics and simulation details. Adjust these parameters as needed.

## Results

Put option, S0 = 100, σ = 20%, r = 5%, T = 1 year, 50 exercise dates. Reference: Bermudan CRR tree with the same exercise dates (5,000 steps).

| Method | Price (K = 100) | Error vs CRR | Over K ∈ [85, 120] | Over σ ∈ [5%, 25%] | Time |
|---|---|---|---|---|---|
| CRR Bermudan (benchmark) | 6.0785 | – | – | – | < 1 s |
| LSM, 100k paths, m = 4 Laguerre | 6.0974 ± 0.0227 | +0.31% | max 0.79% | max 0.44% | 0.4 s |
| Quantization, 200 points/date | 6.0940 | +0.25% | max 2.39% (deep OTM, K = 85) | max 0.26% | ~15 s |
| European (Black-Scholes) | 5.5735 | – | – | – | – |

- **Validation of LSM against Longstaff & Schwartz (2001), Table 1**: within 0.4% of the CRR benchmark on all six cases tested (e.g. S0 = 36, σ = 20%, T = 1: 4.4751 ± 0.0093 vs 4.478 in the paper).
- **Early-exercise premium** captured: 6.08 (Bermudan) vs 5.57 (European) at the money.
- **Quantization** converges to the benchmark as the grid is refined: error +10.8% (20 points), +2.6% (50), +0.8% (100), +0.25% (200), +0.08% (400). The error is largest for deep out-of-the-money puts, where equally likely cells are coarse in the tail.
- **LSM basis size**: m = 1 → −8.2%, m = 2 → −1.2%, m ≥ 3 → within Monte Carlo error.

### Least Squares Monte Carlo vs benchmark

The LSM price matches the Bermudan CRR tree across strikes, and sits above the European price (early-exercise premium):

![LSM vs CRR across strikes](images/lsm_vs_strike.png)

Convergence in the number of paths N (95% confidence intervals):

![LSM convergence in N](images/lsm_convergence.png)

### Quantization tree vs benchmark

![Quantization vs CRR across strikes and volatilities](images/quantization_vs_crr.png)

The relative error decreases steadily as the grid is refined:

![Quantization error vs number of points](images/quantization_convergence.png)

## Corrections

The first version (March 2024, described in `Report.pdf`) under-priced with LSM and gave inconsistent quantization prices. The causes, found and fixed:

1. **r = 0 made the study degenerate.** Without interest rates or dividends, early exercise of a put is never optimal (the payoff is a submartingale), so the American put equals the European put. Prices now use r = 5% with discounting.
2. **LSM regression basis.** `x^k` for k < 30 with x ≈ 100 gives regressors up to 10^58 (overflow, ill-conditioned problem solved with a Powell optimiser). Replaced by m = 4 Laguerre polynomials of S/K, solved exactly by least squares on in-the-money paths.
3. **LSM stopping-time bug.** In `compute_premium`, the loop updating the stopping times was outside the backward loop, and the helper functions read the *global* `alpha` and `tau` from a different simulation. The exercise rule applied to each path came from other paths, which was the main cause of the under-pricing. The exercise decision at the last date before maturity was also never taken.
4. **Quantization transition matrix.** `scipy.integrate.dblquad` passes `(inner, outer)` to the integrand; the variables were swapped and a spurious Jacobian factor was added, so the rows of the transition matrix summed to ~10^-114 instead of 1. There were also index errors (`[i]` instead of `[j]`) and wrong cell bounds. The matrix is now computed from the independence of S_t and S_{t+1}/S_t, with Gauss-Legendre quadrature, and its rows sum to 1.
5. **Benchmark.** The "binomial model" was a *European* tree with only 10 steps. It is now an American/Bermudan CRR tree with thousands of steps.


