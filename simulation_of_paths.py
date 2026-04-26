# SIMULATE PATHS FOR AN APPROXIMATED SOLUTION X_T
import numpy as np
import matplotlib.pyplot as plt

# set up the parameters
num_paths = 15        # Number of paths to simulate
T = 1.0               # Total time to simulate
N = 1000              # Number of time steps
dt = T / N            # Size of each time step
t = np.linspace(0, T, N+1)  # Time vector for plotting
X0 = 10

# Icreate empty arrays to hold the values
X = np.zeros((N+1, num_paths))
X[0] = X0

# simulate the paths using the Milstein method for each time interval $\Delta(T)$
for i in range(N):
    # Generate random noise
    dW = np.random.normal(0, np.sqrt(dt), num_paths)

    # compute the necessary terms in the Milstein scheme
    drift = 1.0 * dt
    diffusion = 2.0 * np.sqrt(np.maximum(X[i], 0)) * dW
    milstein_correction = 1.0 * (dW**2 - dt)

    # update the next step
    X[i+1] = X[i] + drift + diffusion + milstein_correction

# Plotting the chart
plt.figure(figsize=(10, 6))
plt.plot(t, X, lw=1.5)
plt.xlabel('Time (t)', fontsize=12)
plt.title(r'Simulated Paths for the SDE $dX_t = dt + 2\sqrt{X_t}dW_t$')
plt.ylabel('$X_t$', fontsize=12)
plt.grid(True, alpha=0.3)
plt.xlim(0, T)

# Show the plot
plt.show()
