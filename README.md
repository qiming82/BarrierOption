Test implementation of Barrier option with stochastic local model

Used Heston CEV-like model
\begin{eqnarray}
dS_t &= rS_t dt + \sqrt{V_t}\sigma(S_t,t) S_t dW_t^s,\\
dV_t &= \kappa (\theta - V_t) + \xi\sqrt{V_t} dW_t^v
\end{eqnarray}

Smart pointer is used for rng (auto delete to avoid mem leak), seed is provided as input for the replication purpose.

First few functionalities:
- simple monte carlo for price
- greeks including delta, gamma, vega (to be extended if necessary)
- check a range of outcomes when varying intial values (S0, V0)

To do:
- to calculate an implied vol
- to include flag to calculate American option

Output of first few tests:
qimingwang@Qimings-Mac-mini barrierOption_TD % make            
g++ -O3 -Wall -std=c++11 -Iinclude -c src/barrier_analyzer.cpp -o src/barrier_analyzer.o
g++ -O3 -Wall -std=c++11 -Iinclude -o barrier_option main.o src/barrier_analyzer.o src/barrier_option.o
qimingwang@Qimings-Mac-mini barrierOption_TD % ./barrier_option
Up-and-Out Call Option Price: 3.42639
Delta: 0.788809
Gamma: 521.17
Vega: 6.31112
Time cost in secs:      15.8067 seconds
start calculating range ...
computing range ... step: 0
computing range ... step: 1
computing range ... step: 2
computing range ... step: 3
computing range ... step: 4
computing range ... step: 5
computing range ... step: 6
computing range ... step: 7
computing range ... step: 8
computing range ... step: 9
computing range ... step: 10
computing range ... step: 11
computing range ... step: 12
computing range ... step: 13
computing range ... step: 14
computing range ... step: 15
computing range ... step: 16
computing range ... step: 17
computing range ... step: 18
computing range ... step: 19
computing range ... step: 20
computing range ... step: 21
computing range ... step: 22
computing range ... step: 23
computing range ... step: 24
computing range ... step: 25
computing range ... step: 26
computing range ... step: 27
computing range ... step: 28
computing range ... step: 29
computing range ... step: 30
computing range ... step: 31
computing range ... step: 32
computing range ... step: 33
computing range ... step: 34
computing range ... step: 35
computing range ... step: 36
computing range ... step: 37
computing range ... step: 38
computing range ... step: 39
Greeks data saved to greeks.csv
Time cost for range test in secs:      610.087 seconds
