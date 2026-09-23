"""American / Bermudan option pricers under Black-Scholes dynamics.

dS_t = r S_t dt + sigma S_t dW_t  (risk-neutral, no dividends)

- black_scholes : closed-form European price
- crr           : Cox-Ross-Rubinstein binomial tree (European, American or Bermudan)
- lsm           : Longstaff-Schwartz least-squares Monte Carlo
- quantization  : quantization tree (Bally-Pages) with exact transition probabilities
"""
import numpy as np
from numpy.polynomial import laguerre
from scipy.stats import norm


def payoff(S, K, option_type='put'):
    if option_type == 'call':
        return np.maximum(S - K, 0.0)
    if option_type == 'put':
        return np.maximum(K - S, 0.0)
    raise ValueError("option_type must be 'call' or 'put'")


def black_scholes(S, K, T, r, sigma, option_type='put'):
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    if option_type == 'call':
        return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)


def crr(S0, K, T, r, sigma, N, option_type='put', style='american', exercise_dates=None):
    """CRR binomial tree with N steps.

    style: 'european', 'american' (exercise at every step) or 'bermudan'
    (exercise only at the `exercise_dates` equally spaced dates; N must be a multiple of it).
    """
    dt = T / N
    u = np.exp(sigma * np.sqrt(dt))
    d = 1 / u
    p = (np.exp(r * dt) - d) / (u - d)
    disc = np.exp(-r * dt)
    if style == 'bermudan':
        if N % exercise_dates:
            raise ValueError("N must be a multiple of exercise_dates")
        step = N // exercise_dates

    def nodes(i):
        return S0 * u ** np.arange(i, -1, -1) * d ** np.arange(0, i + 1)

    V = payoff(nodes(N), K, option_type)
    for i in range(N - 1, -1, -1):
        V = disc * (p * V[:-1] + (1 - p) * V[1:])
        if style == 'american' or (style == 'bermudan' and i % step == 0):
            V = np.maximum(V, payoff(nodes(i), K, option_type))
    return V[0]


def simulate_paths(S0, T, r, sigma, L, N, seed=None):
    """N paths (N/2 antithetic pairs) of S at the dates 0, T/L, ..., T. Shape (N, L+1)."""
    rng = np.random.default_rng(seed)
    dt = T / L
    G = rng.standard_normal((N // 2, L))
    G = np.vstack([G, -G])
    log_increments = (r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * G
    log_S = np.hstack([np.zeros((G.shape[0], 1)), np.cumsum(log_increments, axis=1)])
    return S0 * np.exp(log_S)


def lsm(S0, K, T, r, sigma, L, N, m=4, option_type='put', seed=None, return_stderr=False):
    """Longstaff-Schwartz price of a Bermudan option exercisable at the L dates T/L, ..., T.

    The continuation value is regressed, on in-the-money paths only, on the first m
    Laguerre polynomials of the normalised spot S/K. Normalising keeps the regressors
    O(1) and the Laguerre family keeps the design matrix well conditioned; the original
    basis x^k (x ~ 100, k up to 29) overflowed and made the regression meaningless.
    """
    S = simulate_paths(S0, T, r, sigma, L, N, seed)
    Z = payoff(S, K, option_type)
    disc = np.exp(-r * T / L)

    cashflow = Z[:, L].copy()  # value of following tau_{j+1}, seen at date j+1
    for j in range(L - 1, 0, -1):
        cashflow *= disc
        itm = Z[:, j] > 0
        if itm.sum() <= m:
            continue
        basis = laguerre.lagvander(S[itm, j] / K, m - 1)
        coef, *_ = np.linalg.lstsq(basis, cashflow[itm], rcond=None)
        exercise = Z[itm, j] > basis @ coef
        idx = np.flatnonzero(itm)[exercise]
        cashflow[idx] = Z[idx, j]
    cashflow *= disc

    price = max(Z[0, 0], cashflow.mean())
    if return_stderr:
        return price, cashflow.std(ddof=1) / np.sqrt(N)
    return price


def _strata(n):
    """Equal-probability strata of N(0,1): boundaries (n+1,) and conditional means (n,)."""
    b = norm.ppf(np.linspace(0, 1, n + 1))
    return b, (norm.pdf(b[:-1]) - norm.pdf(b[1:])) * n


def quantization_grid(S0, t, r, sigma, n):
    """Quantize S_t on n equally likely cells.

    Returns cell boundaries in S (n+1,) and the quantizer of each cell, chosen as the
    conditional mean E[S_t | S_t in cell] so that the quantized process keeps E[S_t].
    """
    b, _ = _strata(n)
    mu = np.log(S0) + (r - 0.5 * sigma**2) * t
    s = sigma * np.sqrt(t)
    bounds = np.exp(mu + s * b)
    # E[e^{mu+sG} 1{a<G<b}] = e^{mu+s^2/2} (Phi(b-s) - Phi(a-s))
    points = np.exp(mu + 0.5 * s**2) * np.diff(norm.cdf(b - s)) * n
    return bounds, points


def transition_matrix(S0, t, dt, r, sigma, n, n_quad=40):
    """pi[i, j] = P(S_{t+dt} in cell j | S_t in cell i), for equal-probability grids.

    S_{t+dt} = S_t * R with R lognormal and independent of S_t, so
    pi[i, j] = n * int_{cell i} f_{S_t}(y) P(R in [a_j / y, b_j / y]) dy.
    The outer integral is done by Gauss-Legendre in the Gaussian variable of S_t.
    """
    b, _ = _strata(n)
    mu_t = np.log(S0) + (r - 0.5 * sigma**2) * t
    s_t = sigma * np.sqrt(t)
    next_bounds, _ = quantization_grid(S0, t + dt, r, sigma, n)
    m_R = (r - 0.5 * sigma**2) * dt
    s_R = sigma * np.sqrt(dt)

    x, w = np.polynomial.legendre.leggauss(n_quad)
    # finite truncation of the two infinite outer strata (mass beyond 8 sd is < 1e-15)
    lo = np.maximum(b[:-1], -8.0)
    hi = np.minimum(b[1:], 8.0)
    g = 0.5 * (hi - lo)[:, None] * x[None, :] + 0.5 * (hi + lo)[:, None]  # (n, n_quad)
    weights = 0.5 * (hi - lo)[:, None] * w[None, :] * norm.pdf(g)           # (n, n_quad)
    log_y = mu_t + s_t * g

    with np.errstate(divide='ignore'):
        log_next = np.log(next_bounds)                                      # -inf, ..., +inf
    cdf = norm.cdf((log_next[None, None, :] - log_y[:, :, None] - m_R) / s_R)
    cell_prob = np.diff(cdf, axis=2)                                        # (n, n_quad, n)
    pi = n * np.einsum('iq,iqj->ij', weights, cell_prob)
    return pi / pi.sum(axis=1, keepdims=True)


def quantization(S0, K, T, r, sigma, L, n=100, option_type='put'):
    """Quantization tree price of a Bermudan option exercisable at the L dates T/L, ..., T."""
    dt = T / L
    disc = np.exp(-r * dt)

    _, points = quantization_grid(S0, T, r, sigma, n)
    V = payoff(points, K, option_type)
    for k in range(L - 1, 0, -1):
        _, points = quantization_grid(S0, k * dt, r, sigma, n)
        pi = transition_matrix(S0, k * dt, dt, r, sigma, n)
        V = np.maximum(payoff(points, K, option_type), disc * pi @ V)
    # from t=0 (deterministic S0) every cell of date 1 has probability 1/n
    return max(payoff(S0, K, option_type), disc * V.mean())
