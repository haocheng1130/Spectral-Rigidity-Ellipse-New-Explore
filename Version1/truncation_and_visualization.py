# For visualization1,2
import numpy as np
import pickle
import matplotlib.pyplot as plt
from numpy.linalg import eig

# For visualization3
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
print(animation.writers.list())

# Load the precomputed T_D_inverse matrix
T_D_inverse = np.loadtxt("matrix_T_D_inverse_10000x10000.csv", delimiter=',')

print(f"\nT_D_inverse (first 5 rows):")
print(T_D_inverse[:5, :5])

# Load the reduced matrices from the pickle file
with open("reduced_matrices.pkl", "rb") as file:
    reduced_matrices = pickle.load(file)
    print("Reduced matrices loaded.")




# ------------------------------------------------------
# Visualization1: This is to draw all the eigenvalue for all truncations on a single plot.
# ------------------------------------------------------


# Parameters for truncation
min_N = 10
max_N = 5000
step = 10
Ns = list(range(min_N, max_N + 1, step))

# Function to compute spectra for submatrices
def compute_spectrum(matrix, N):
    """
    Truncates the matrix to an N x N submatrix and computes its spectrum (eigenvalues).
    """
    truncated_matrix = matrix[:N, :N]  # Truncate to N x N
    eigenvalues = eig(truncated_matrix)[0]  # Compute eigenvalues
    return eigenvalues

# Iterate over eccentricities and generate plots
for e, matrix in reduced_matrices.items():
    plt.figure(figsize=(12, 8))
    colors = plt.cm.viridis(np.linspace(0, 1, len(Ns)))  # Color map for N steps

    # Compute the matrix required for truncation
    required_matrix = np.dot(T_D_inverse[:max_N, :max_N], matrix)

    for i, N in enumerate(Ns):
        # Compute the spectrum for the N x N submatrix
        spectrum = compute_spectrum(required_matrix, N)
        
        # Plot the spectrum on the complex plane
        plt.scatter(
            np.real(spectrum),  # Real part of eigenvalues (x-axis)
            np.imag(spectrum),  # Imaginary part of eigenvalues (y-axis)
            color=colors[i],
            label=f"N={N}" if i == 0 else "",
            s=10
        )

    # Configure the plot for this eccentricity
    plt.title(f"Eigenvalues in the Complex Plane for Eccentricity e={e}")
    plt.xlabel("Real Part of Eigenvalues")
    plt.ylabel("Imaginary Part of Eigenvalues")
    plt.colorbar(plt.cm.ScalarMappable(cmap="viridis"), label="Truncation Size (N)")
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)  # Horizontal line
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)  # Vertical line
    plt.grid(True)
    plt.savefig(f"spectra_complex_plane_e_{e:.2f}.png", dpi=300)  # Save the plot
    plt.show()


# ------------------------------------------------------
# Visualization2: For each truncation, we generate a plot.
# ------------------------------------------------------

# Parameters for truncation
min_N = 10
max_N = 5000
step = 10
Ns = list(range(min_N, max_N + 1, step))

# Function to compute spectra for submatrices
def compute_spectrum(matrix, N):
    """
    Truncates the matrix to an N x N submatrix and computes its spectrum (eigenvalues).
    """
    truncated_matrix = matrix[:N, :N]  # Truncate to N x N
    eigenvalues = eig(truncated_matrix)[0]  # Compute eigenvalues
    return eigenvalues

# Iterate over eccentricities and generate plots for each truncation size
for e, matrix in reduced_matrices.items():
    print(f"Processing eccentricity e={e}")

    # Compute the matrix required for truncation
    required_matrix = np.dot(T_D_inverse[:max_N, :max_N], matrix)

    for N in Ns:
        # Compute the spectrum for the N x N submatrix
        spectrum = compute_spectrum(required_matrix, N)

        # Create a plot for the complex plane
        plt.figure(figsize=(10, 8))
        plt.scatter(
            np.real(spectrum),  # Real part (x-axis)
            np.imag(spectrum),  # Imaginary part (y-axis)
            c='blue',
            s=20,
            label=f"Eigenvalues for N={N}"
        )

        # Configure the plot
        plt.title(f"Eigenvalues in the Complex Plane\nEccentricity e={e}, Truncation Size N={N}")
        plt.xlabel("Real Part of Eigenvalues")
        plt.ylabel("Imaginary Part of Eigenvalues")
        plt.axhline(0, color='gray', linestyle='--', linewidth=0.5)  # Horizontal line
        plt.axvline(0, color='gray', linestyle='--', linewidth=0.5)  # Vertical line
        plt.grid(True)
        plt.legend(loc='upper right')

        # Save the plot for this truncation size
        filename = f"spectra_complex_plane_e_{e:.2f}_N_{N}.png"
        plt.savefig(filename, dpi=300)
        plt.close()  # Close the figure to save memory
        print(f"Saved plot for N={N}, eccentricity e={e} to {filename}")



# ------------------------------------------------------
# Visualization3: Animation
# ------------------------------------------------------

# Parameters for truncation
min_N = 10
max_N = 5000
step = 10
Ns = list(range(min_N, max_N + 1, step))

# Function to compute spectra for submatrices
def compute_spectrum(matrix, N):
    """
    Truncates the matrix to an N x N submatrix and computes its spectrum (eigenvalues).
    """
    truncated_matrix = matrix[:N, :N]  # Truncate to N x N
    eigenvalues = eig(truncated_matrix)[0]  # Compute eigenvalues
    return eigenvalues

# Select the eccentricity and matrix for animation
e = list(reduced_matrices.keys())[0]  # Example: Use the first eccentricity
matrix = reduced_matrices[e]
print(f"Animating for eccentricity e={e}")

# Compute the matrix required for truncation
required_matrix = np.dot(T_D_inverse[:max_N, :max_N], matrix)

# Prepare the figure
fig, ax = plt.subplots(figsize=(10, 8))
scatter = ax.scatter([], [], c='blue', s=20)
ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)
ax.set_title(f"Eigenvalues in the Complex Plane\nEccentricity e={e}")
ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)
ax.set_xlabel("Real Part of Eigenvalues")
ax.set_ylabel("Imaginary Part of Eigenvalues")
ax.grid(True)

# Initialize the plot
def init():
    scatter.set_offsets(np.array([[0, 0]]))  # Initialize with a dummy point at (0, 0)
    return scatter,

# Update function for animation
def update(frame):
    N = Ns[frame]
    spectrum = compute_spectrum(required_matrix, N)
    scatter.set_offsets(np.c_[np.real(spectrum), np.imag(spectrum)])
    ax.set_title(f"Eigenvalues in the Complex Plane\nEccentricity e={e}, Truncation Size N={N}")
    return scatter,

# Create the animation
anim = FuncAnimation(fig, update, frames=len(Ns), init_func=init, blit=False)

# Save as a video
anim.save("eigenvalues_complex_plane2.mp4", writer='ffmpeg', fps=2, dpi=300)
print("Animation saved as eigenvalues_complex_plane2.mp4")
