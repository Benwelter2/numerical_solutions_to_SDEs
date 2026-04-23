''''
This code computes the strong order of convergence for the Euler scheme.
Therefore, we generate 10'000 independent sample paths of a stochastic process,
as well as 10'000 approximated sample paths using the Euler scheme. We then compute
the average error between the approximated and the exact solution. Next, we plot
log_2(N) against -log_2(epsilon) and for a specific N (like 2**4), we found an average error ϵ
and we take the "log"-linear fit of these averaged erros for each N. The slope of this
line equals the order of strong convergence.
'''

# import libraries
import numpy as np
import matplotlib.pyplot as plt

# setting up parameters
T = 1.0
N = np.array([2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**10])
num_paths = 10**4
mu = 5
sigma = 0.6
X_0 = 10.0

# empty array to store the mean errors
epsilon = np.zeros(len(N))

# Loop over grid sizes
for n, N_n in enumerate(N):
    dt = T / N_n

    # initialize Y (the approximated solution) and W_t (to store the BM-value) for each path
    Y = np.full(num_paths, X_0)
    W_T = np.zeros(num_paths) 

    # Loop over time steps
    for i in range(N_n):
        # Generate dW for all the paths at once
        dW = np.sqrt(dt) * np.random.randn(num_paths)
        W_T += dW

        # update the Euler scheme on each time increment
        Y += mu * Y * dt + sigma * Y * dW

    # compute the exact solution of the SDE at time T
    X_T = X_0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)

    # compute the strong error
    error = np.abs(X_T - Y)
    epsilon[n] = np.mean(error)

# compute log_2(N) and -log_2(epsilon) 
log2_N = np.log2(N)
neg_log2_eps = -np.log2(epsilon)

# compute the linear fit line of the averaged erros
line_fit = np.polyfit(log2_N, neg_log2_eps, 1)
print(f"Line fit (Slope, Intercept): {line_fit}")

# Plotting log_2(N) against -log_2(epsilon) with the linear fit line
plt.plot(log2_N, neg_log2_eps, '+', markersize=8, label='Simulation Errors')
plt.plot(log2_N, line_fit[0] * log2_N + line_fit[1], ':', label=f'Linear Fit (Slope: {line_fit[0]:.4f})')
plt.xlabel('$log_2(N)$')
plt.ylabel('$-log_2(\epsilon)$')
plt.title('Strong Convergence of the Euler Scheme')
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0, 11)  
plt.ylim(0, 10)
plt.legend(loc='upper left')
plt.show()
