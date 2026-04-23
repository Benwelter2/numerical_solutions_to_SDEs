''''
This code computes the strong order of convergence for the Milstein scheme.
Therefore, we generate 10'000 independent sample paths of a stochastic process,
as well as 10'000 approximated sample paths using the Milstein scheme. We then compute the absolute
difference between the average of the functional over all simulated paths. Next, we plot
log_2(N) against -log_2(epsilon) for various step sizes. By applying a linear fit to these
log-transformed errors, the resulting slope represents the order of weak convergence.
'''

# import libraries
import numpy as np
import matplotlib.pyplot as plt

# setting up parameters
T = 1.0
N = np.array([2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**10])
J = 10**5 # 
mu = 0.5
sigma = 0.6
X_0 = 10.0

# empty array to store the mean errors
epsilon = np.zeros(len(N))

# Loop over grid sizes
for n, N_n in enumerate(N):
    dt = T / N_n

    # initialize Y (the approximated solution) and W_t (to store the BM-value) for each path
    Y = np.full(J, X_0)
    W_T = np.zeros(J) 

    # Loop over time steps
    for i in range(N_n):
        # Generate dW for all the paths at once
        dW = np.sqrt(dt) * np.random.randn(J)
        W_T += dW

        # update the Milstein scheme on each time increment
        Y += mu * Y * dt + sigma * Y * dW + 0.5 * sigma**2 * Y * (dW**2 - dt)

    # compute the exact solution of the SDE at time T
    X_T = X_0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)

    # compute weak error |E[X_T] - E[Y_T]| using g(x) = x
    epsilon[n] = np.abs(np.mean(X_T) - np.mean(Y))

# compute log_2(N) and -log_2(epsilon)
log2_N = np.log2(N)
neg_log2_eps = -np.log2(epsilon)

# compute the linear fit line
line_fit = np.polyfit(log2_N, neg_log2_eps, 1)
print(f"Line fit (Slope, Intercept): {line_fit}")

# Plotting log_2(N) against -log_2(epsilon) with the linear fit line
plt.plot(log2_N, neg_log2_eps, '+', markersize=8, label='Simulation Errors')
plt.plot(log2_N, line_fit[0] * log2_N + line_fit[1], ':', label=f'Linear Fit (Slope: {line_fit[0]:.4f})')
plt.xlabel('$log_2(N)$')
plt.ylabel('$-log_2(\epsilon)$')
plt.title('Weak Convergence of the Milstein Scheme')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.xlim(0, 11)  
plt.ylim(0, 10)
plt.show()