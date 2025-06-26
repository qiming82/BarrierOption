import numpy as np
import matplotlib.pyplot as plt

# Heston parameters
r = 0.05
kappa = 2.0
theta = 0.04
sigma = 0.3
rho = -0.7
v0 = 0.04

# Option parameters
K = 100
B = 120
T = 1.0

# Grid parameters
Smax, Smin, M = 1.2 * B, 0, 80
vmax, vmin, N = 0.5, 0.0, 40
L = 100  # time steps

# Adaptive S grid: denser around K
def make_adaptive_S_grid(Smin, Smax, M, K, density=5.0):
    # density: higher = more points near K
    x = np.linspace(-1, 1, M)
    S = K + (Smax-K) * np.sinh(density * x) / np.sinh(density)
    S = np.clip(S, Smin, Smax)
    return S

S = make_adaptive_S_grid(Smin, Smax, M, K, density=5.0)
v = np.linspace(vmin, vmax, N)
dt = T / L
dv = v[1] - v[0]
t_grid = np.linspace(0, T, L+1)

# Example Dupire local volatility function (smile)
def sigma_loc(S, t):
    return 0.15 + 0.2 * np.exp(-((S-K)/30)**2)

# Payoff: up-and-out call
def payoff(S):
    return np.maximum(S - K, 0) * (S < B)

# Initial condition
V = np.zeros((M, N))
for i in range(M):
    V[i, :] = payoff(S[i])

# Thomas algorithm for tridiagonal systems
def thomas(a, b, c, d):
    n = len(d)
    c_ = np.zeros(n-1)
    d_ = np.zeros(n)
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    for i in range(1, n-1):
        c_[i] = c[i] / (b[i] - a[i-1] * c_[i-1])
    for i in range(1, n):
        d_[i] = (d[i] - a[i-1] * d_[i-1]) / (b[i] - a[i-1] * c_[i-1])
    x = np.zeros(n)
    x[-1] = d_[-1]
    for i in range(n-2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i+1]
    return x

# ADI time stepping with non-uniform S grid
for l in range(L):
    t = t_grid[L-l]  # backward in time
    # S-direction implicit
    V_temp = np.copy(V)
    for j in range(1, N-1):
        A = np.zeros(M)
        B_ = np.zeros(M)
        C = np.zeros(M)
        D = np.zeros(M)
        for i in range(1, M-1):
            dS_plus = S[i+1] - S[i]
            dS_minus = S[i] - S[i-1]
            dS_avg = 0.5 * (dS_plus + dS_minus)
            sigloc = sigma_loc(S[i], t)
            vfac = v[j] * sigloc**2
            # Second derivative (central, non-uniform)
            a = vfac * S[i]**2 / (dS_minus * dS_avg)
            c = vfac * S[i]**2 / (dS_plus * dS_avg)
            # First derivative (central, non-uniform)
            b = r * S[i] * (1/dS_plus - 1/dS_minus) / 2
            # Tridiagonal coefficients
            A[i] = -dt * (a - b)
            B_[i] = 1 + dt * (a + c + r)
            C[i] = -dt * (c + b)
            D[i] = V[i, j]
        # Boundary conditions
        B_[0] = 1
        C[0] = 0
        D[0] = 0
        A[-1] = 0
        B_[-1] = 1
        D[-1] = 0
        V_temp[:, j] = thomas(A[1:], B_, C[:-1], D)
    # v-direction implicit
    for i in range(1, M-1):
        A = np.zeros(N)
        B_ = np.zeros(N)
        C = np.zeros(N)
        D = np.zeros(N)
        for j in range(1, N-1):
            c = 0.5 * sigma**2 * v[j] / dv**2
            d = kappa * (theta - v[j]) / (2 * dv)
            A[j] = -dt * (c - d)
            B_[j] = 1 + dt * (2 * c + r)
            C[j] = -dt * (c + d)
            D[j] = V_temp[i, j]
        # Boundary conditions
        B_[0] = 1
        C[0] = 0
        D[0] = 0
        A[-1] = 0
        B_[-1] = 1
        D[-1] = 0
        V[i, :] = thomas(A[1:], B_, C[:-1], D)
    # Barrier condition
    V[S >= B, :] = 0

# Extract price and Greeks at v0
j0 = np.argmin(np.abs(v - v0))
price = V[:, j0]
# Greeks on non-uniform grid
delta = np.zeros(M)
gamma = np.zeros(M)
for i in range(1, M-1):
    delta[i] = (price[i+1] - price[i-1]) / (S[i+1] - S[i-1])
    gamma[i] = 2 * ((price[i+1] - price[i]) / (S[i+1] - S[i]) - (price[i] - price[i-1]) / (S[i] - S[i-1])) / (S[i+1] - S[i-1])
delta[0] = delta[1]
delta[-1] = delta[-2]
gamma[0] = gamma[1]
gamma[-1] = gamma[-2]
vega = (V[:, j0+1] - V[:, j0-1]) / (v[j0+1] - v[j0-1])

# Plotting
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
axs[0, 0].plot(S, price, label='Option Price')
axs[0, 0].axvline(K, color='r', linestyle='--', label='Strike (K=100)')
axs[0, 0].axvline(B, color='g', linestyle='--', label='Barrier (B=120)')
axs[0, 0].set_ylabel('Price')
axs[0, 0].legend()

axs[0, 1].plot(S, delta, label='Delta')
axs[0, 1].axvline(K, color='r', linestyle='--')
axs[0, 1].axvline(B, color='g', linestyle='--')
axs[0, 1].set_ylabel('Delta')
axs[0, 1].legend()

axs[1, 0].plot(S, gamma, label='Gamma')
axs[1, 0].axvline(K, color='r', linestyle='--')
axs[1, 0].axvline(B, color='g', linestyle='--')
axs[1, 0].set_ylabel('Gamma')
axs[1, 0].legend()

axs[1, 1].plot(S, vega, label='Vega')
axs[1, 1].axvline(K, color='r', linestyle='--')
axs[1, 1].axvline(B, color='g', linestyle='--')
axs[1, 1].set_ylabel('Vega')
axs[1, 1].set_xlabel('Stock Price (S)')
axs[1, 1].legend()

fig.suptitle('Greeks for Up-and-Out Call Option (Heston + Dupire Local Vol, Adaptive Grid)')
plt.tight_layout()
plt.show()