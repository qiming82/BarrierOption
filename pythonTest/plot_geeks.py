import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('greeks.csv')

fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharex=True)

axs[0,0].plot(data['StockPrice'], data['Price'], label='Option Price')
axs[0,0].axvline(x=100, color='r', linestyle='--', label='Strike (K=100)')
axs[0,0].axvline(x=120, color='g', linestyle='--', label='Barrier (B=120)')
axs[0,0].set_ylabel('Price')
axs[0,0].legend()
axs[0,0].grid(True)

axs[0,1].plot(data['StockPrice'], data['Delta'], label='Delta')
axs[0,1].axvline(x=100, color='r', linestyle='--')
axs[0,1].axvline(x=120, color='g', linestyle='--')
axs[0,1].set_ylabel('Delta')
axs[0,1].legend()
axs[0,1].grid(True)

axs[1,0].plot(data['StockPrice'], data['Gamma'], label='Gamma')
axs[1,0].axvline(x=100, color='r', linestyle='--')
axs[1,0].axvline(x=120, color='g', linestyle='--')
axs[1,0].set_ylabel('Gamma')
axs[1,0].legend()
axs[1,0].grid(True)

axs[1,1].plot(data['StockPrice'], data['Vega'], label='Vega')
axs[1,1].axvline(x=100, color='r', linestyle='--')
axs[1,1].axvline(x=120, color='g', linestyle='--')
axs[1,1].set_xlabel('Stock Price (S)')
axs[1,1].set_ylabel('Vega')
axs[1,1].legend()
axs[1,1].grid(True)

plt.suptitle('Greeks for Up-and-Out Call Option (SLV Model)')
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('greeks_plot.png')
plt.show()