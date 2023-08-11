import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# User inputs
L = float(input("Enter the length of the region (L): "))
N = int(input("Enter the number of spatial points (N): "))
T = float(input("Enter the total time (T): "))
num_steps = int(input("Enter the number of time steps: "))
sigma = float(input("Enter the width of the wave packet (sigma): "))
x0 = float(input("Enter the initial position of the wave packet (x0): "))

# Discretization
dx = L / (N - 1)
dt = T / num_steps

# Initialize the wave function
x = np.linspace(0, L, N)
psi = np.exp(-0.5 * ((x - x0) / sigma) ** 2 + 1j * x * 5)  # Example initial wave packet


# Define the potential energy function (harmonic oscillator potential in this case)
def potential(x):
    return 0.5 * x ** 2


# Create the sparse matrix for the Hamiltonian
diagonal = np.ones(N) * (1 / dx ** 2) + potential(x)
off_diagonal = -0.5 * np.ones(N - 1) * (1 / dx ** 2)
H = diags([off_diagonal, diagonal, off_diagonal], [-1, 0, 1], shape=(N, N))

# Time evolution using Crank-Nicolson method
for step in range(num_steps):
    psi = np.linalg.solve(np.eye(N) - 0.5j * dt * H, psi)

    # Plot the probability density every few steps
    if step % (num_steps // 10) == 0:
        plt.plot(x, np.abs(psi) ** 2, label=f"Step {step}")
        plt.xlabel('Position')
        plt.ylabel('Probability Density')
        plt.legend()
        plt.title(f'Time evolution at step {step}')
        plt.show()

# Plot the final probability density
plt.plot(x, np.abs(psi) ** 2, label=f"Final Step")
plt.xlabel('Position')
plt.ylabel('Probability Density')
plt.legend()
plt.title(f'Final Time evolution')
plt.show()
