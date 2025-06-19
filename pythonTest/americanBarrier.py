import numpy as np
import math
from collections import defaultdict
import time

class BarrierOption:
    def __init__(self, config_file):
        # Read config file
        params = defaultdict(lambda: None)
        with open(config_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    key, _, value = line.strip().partition('=')
                    params[key.strip()] = value.strip()

        # Parse parameters
        try:
            self.S0 = float(params['S0'])
            self.K = float(params['K'])
            self.B = float(params['B'])
            self.T = float(params['T'])
            self.r = float(params['r'])
            self.v0 = float(params['v0'])
            self.kappa = float(params['kappa'])
            self.theta = float(params['theta'])
            self.xi = float(params['xi'])
            self.beta = float(params['beta'])
            self.rho = float(params['rho'])
            self.num_simulations = int(params['numSimulations'])
            self.num_steps = int(params['numSteps'])
            self.seed = int(params['seed']) if params['seed'] else np.random.randint(0, 2**32-1)
            self.option_type = params['optionType']
            self.barrier_type = params['barrierType']
            self.exercise_style = params['exerciseStyle']
        except (KeyError, ValueError) as e:
            raise RuntimeError(f"Invalid config: {e}")

        # Validate parameters
        if any(x <= 0 for x in [self.S0, self.K, self.B, self.T, self.kappa, self.theta, self.xi]):
            raise ValueError("Parameters must be positive (except rho)")
        if self.v0 < 0:
            raise ValueError("Initial variance must be non-negative")
        if abs(self.rho) > 1:
            raise ValueError("Correlation must be between -1 and 1")
        if self.barrier_type == 'UpAndOut' and self.B <= self.S0:
            raise ValueError("Up-and-out barrier must be above initial price")
        if self.barrier_type == 'DownAndOut' and self.B >= self.S0:
            raise ValueError("Down-and-out barrier must be below initial price")
        if self.option_type not in ['Call', 'Put']:
            raise ValueError("Invalid optionType")
        if self.barrier_type not in ['UpAndOut', 'DownAndOut']:
            raise ValueError("Invalid barrierType")
        if self.exercise_style not in ['European', 'American']:
            raise ValueError("Invalid exerciseStyle")

        # Initialize RNG
        self.rng = np.random.default_rng(self.seed)

        # Define model functions
        self.mu_v = lambda v, t: self.kappa * (self.theta - v)
        self.sigma_v = lambda v, t: self.xi * math.sqrt(max(v, 0.0))
        self.local_vol = lambda S, t: 0.2 * (S / 100.0) ** self.beta

    def generate_correlated_normals(self):
        """Generate two correlated standard normal variables."""
        z1 = self.rng.standard_normal()
        z2 = self.rho * z1 + math.sqrt(1.0 - self.rho ** 2) * self.rng.standard_normal()
        return z1, z2

    def simulate_path(self):
        """Simulate path for European option, return payoff."""
        dt = self.T / self.num_steps
        S = self.S0
        v = self.v0
        knocked_out = False

        for i in range(self.num_steps):
            t = i * dt
            z_s, z_v = self.generate_correlated_normals()
            
            v += self.mu_v(v, t) * dt + self.sigma_v(v, t) * math.sqrt(dt) * z_v
            v = max(v, 0.0)
            
            vol = math.sqrt(max(v, 0.0) * self.local_vol(S, t) ** 2)
            S *= math.exp((self.r - 0.5 * vol ** 2) * dt + vol * math.sqrt(dt) * z_s)

            if self.barrier_type == 'UpAndOut' and S >= self.B:
                knocked_out = True
                break
            elif self.barrier_type == 'DownAndOut' and S <= self.B:
                knocked_out = True
                break

        payoff = 0.0
        if not knocked_out:
            if self.option_type == 'Call':
                payoff = max(S - self.K, 0.0)
            else:
                payoff = max(self.K - S, 0.0)
        return payoff

    def simulate_path_with_storage(self):
        """Simulate path for American option, return payoff and path."""
        dt = self.T / self.num_steps
        S = self.S0
        v = self.v0
        knocked_out = False
        path = np.zeros((self.num_steps + 1, 2))
        path[0] = [S, v]

        for i in range(self.num_steps):
            t = i * dt
            z_s, z_v = self.generate_correlated_normals()
            
            v += self.mu_v(v, t) * dt + self.sigma_v(v, t) * math.sqrt(dt) * z_v
            v = max(v, 0.0)
            
            vol = math.sqrt(max(v, 0.0) * self.local_vol(S, t) ** 2)
            S *= math.exp((self.r - 0.5 * vol ** 2) * dt + vol * math.sqrt(dt) * z_s)

            if self.barrier_type == 'UpAndOut' and S >= self.B:
                knocked_out = True
                break
            elif self.barrier_type == 'DownAndOut' and S <= self.B:
                knocked_out = True
                break

            path[i + 1] = [S, v]

        payoff = 0.0
        if not knocked_out:
            if self.option_type == 'Call':
                payoff = max(S - self.K, 0.0)
            else:
                payoff = max(self.K - S, 0.0)
        return payoff, path

    def price_european(self):
        """Calculate price for European option using Monte Carlo."""
        payoffs = np.array([self.simulate_path() for _ in range(self.num_simulations)])
        return math.exp(-self.r * self.T) * np.mean(payoffs)

    def price_american(self):
        """Calculate price for American option using Longstaff-Schwartz."""
        dt = self.T / self.num_steps
        paths = np.zeros((self.num_simulations, self.num_steps + 1, 2))
        cash_flows = np.zeros(self.num_simulations)

        # Simulate paths
        for i in range(self.num_simulations):
            payoff, path = self.simulate_path_with_storage()
            paths[i] = path
            cash_flows[i] = payoff

        # Backward induction
        for t in range(self.num_steps - 1, 0, -1):
            X = []
            Y = []
            indices = []
            for i in range(self.num_simulations):
                S = paths[i, t, 0]
                if self.barrier_type == 'UpAndOut' and S >= self.B:
                    continue
                if self.barrier_type == 'DownAndOut' and S <= self.B:
                    continue
                exercise = max(S - self.K, 0.0) if self.option_type == 'Call' else max(self.K - S, 0.0)
                if exercise > 0:
                    X.append(S)
                    Y.append(cash_flows[i] * math.exp(-self.r * dt * (self.num_steps - t)))
                    indices.append(i)

            if not X:
                continue

            # Quadratic regression: E[Y] = a + b*S + c*S^2
            X = np.array(X)
            Y = np.array(Y)
            A = np.vstack([np.ones_like(X), X, X**2]).T
            try:
                coeffs, _, _, _ = np.linalg.lstsq(A, Y, rcond=None)
            except np.linalg.LinAlgError:
                continue

            # Update cash flows
            for idx in indices:
                S = paths[idx, t, 0]
                exercise = max(S - self.K, 0.0) if self.option_type == 'Call' else max(self.K - S, 0.0)
                continuation = coeffs[0] + coeffs[1] * S + coeffs[2] * S**2 if len(coeffs) == 3 else 0.0
                if exercise > continuation:
                    cash_flows[idx] = exercise * math.exp(self.r * dt * (self.num_steps - t))

        return np.mean(cash_flows)

    def price(self):
        """Calculate option price based on exercise style."""
        if self.exercise_style == 'American':
            return self.price_american()
        return self.price_european()

if __name__ == '__main__':
    try:
        t1 = time.time()
        option = BarrierOption('config.txt')
        price = option.price()
        print(f"{option.barrier_type} {option.option_type} Option Price: {price:.4f}")
        dt = time.time() - t1
        print("elapsed time:", dt)
    except Exception as e:
        print(f"Error: {e}")