# SIMULATE POSSIBLE PATHS FOR A SOLUTION X_T
import numpy as np
import matplotlib.pyplot as plt

# 1. Set up the parameters
num_paths = 15   # Number of paths to simulate
T = 1.0               # Total time to simulate
N = 1000              # Number of time steps
dt = T / N            # Size of each time step
t = np.linspace(0, T, N+1)  # Time vector for plotting

# Initial condition
X0 = 10

# 2. Initialize the array to hold our paths
X = np.zeros((N+1, num_paths))
X[0] = X0

# 3. Simulate the paths using the Milstein method
for i in range(N):
    # Generate random standard normal noise
    dW = np.random.normal(0, np.sqrt(dt), num_paths)

    # Calculate standard Euler-Maruyama terms
    drift = 1.0 * dt
    diffusion = 2.0 * np.sqrt(np.maximum(X[i], 0)) * dW

    # Calculate the Milstein correction term
    # 1/2 * b(X) * b'(X) simplifies to exactly 1.0 for this specific SDE
    milstein_correction = 1.0 * (dW**2 - dt)

    # Update the next step
    X[i+1] = X[i] + drift + diffusion + milstein_correction

# 4. Plot the "Spaghetti Plot"
plt.figure(figsize=(10, 6))
plt.plot(t, X, lw=1.5)

# Formatting the chart
# Note the "r" before the string to correctly render the LaTeX in matplotlib
plt.xlabel('Time (t)', fontsize=12)
plt.title(r'Simulated Paths for the SDE $dX_t = dt + 2\sqrt{X_t}dW_t$')
plt.ylabel('$X_t$', fontsize=12)
plt.grid(True, alpha=0.3)
plt.xlim(0, T)

# Show the plot
plt.show()
