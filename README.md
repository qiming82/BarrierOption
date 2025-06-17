Test implementation of Barrier option with stochastic local model

Used Heston CEV-like model
\begin{eqnarray}
dS_t &= rS_t dt + \sqrt{V_t}\sigma(S_t,t) S_t dW_t^s,\\
dV_t &= \kappa (\theta - V_t) + \xi\sqrt{V_t} dW_t^v
\end{eqnarray}

\begin{equation}
    \sigma = 0.2 \left(\frac{S}{S_0}\right)^\beta
\end{equation}
if \sigma = 1, it is just Heston model.


Smart pointer is used for rng (auto delete to avoid mem leak), seed is provided as input for the replication purpose.

First few functionalities:
- simple monte carlo for price
- greeks including delta, gamma, vega (to be extended if necessary)
- check a range of outcomes when varying intial values (S0, V0)
- to add convergence test

To do:
- to include code to calculate American option?
