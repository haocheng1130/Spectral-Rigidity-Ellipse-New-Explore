from mpmath import mp, mpf, sin, pi, nstr
import pandas as pd
import pickle
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------------------------
# Part 0: Helper functions
# ------------------------------------------------------------------------------------------------------------

def mat_mult(A, B):
    """
    Multiply two square matrices A and B of dimension n x n.
    Returns the resulting matrix as a list of lists.
    """
    n = len(A)
    # Initialize the result matrix with zeros (using high-precision mpf)
    C = [[mp.mpf(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = mp.mpf(0)
            for k in range(n):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C


# ------------------------------------------------------------------------------------------------------------
# Part 1: Load all the matrices we need to consider
# ------------------------------------------------------------------------------------------------------------

with open("matrix_T_D_inverse.pkl", "rb") as f:
    T_D_inverse = pickle.load(f)  

with open("T_matrices.pkl", "rb") as f:
    T_matrices = pickle.load(f)  


# ------------------------------------------------------------------------------------------------------------
# Part 2: Plot all the eigenvalues on the complex plane
# ------------------------------------------------------------------------------------------------------------

list_e = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
# We only consider 10 eccentricities for simplicity
# Running this code will takes a very long time, so I suggest to run the code for one eccentricity at a time
# If you want to run the code for other e, please check T_matrices has that key


for e in list_e:
    print(f"Processing eccentricity e={e}")
    T = T_matrices[e]
    test_matrix = mat_mult(T_D_inverse, T)

    for size in range(50, 300, 10):
        # Extract top-left 'size' x 'size' submatrix
        submatrix = mp.matrix([row[:size] for row in test_matrix[:size]])
        
        # Compute eigenvalues (for a real symmetric or complex Hermitian matrix)
        # If your matrix is definitely Hermitian, you can set eigvals_only=False
        # to also get eigenvectors. Here we only need the eigenvalues.
        eigenvalues, _ = mp.eig(submatrix)
        
        # Separate real and imaginary parts for plotting
        real_parts = [ev.real for ev in eigenvalues]
        imag_parts = [ev.imag for ev in eigenvalues]
        
        # Plot the eigenvalues on the complex plane
        # We used high precision for other stuff, if you only want these pictures, please use lower precision
        plt.figure()
        plt.scatter(real_parts, imag_parts, c='blue', marker='o', s=10)
        plt.title(f"Eigenvalues in the Complex Plane\nEccentricity e={e}, Truncation Size={size}")
        plt.xlabel("Real Part")
        plt.ylabel("Imaginary Part")
        plt.grid(True)
        
        # Save the figure to a PNG file
        plt.savefig(f"spectrum_e_{e:.2f}_size_{size}.png")
        plt.close()
        
        print(f"Saved spectrum_e_{e:.2f}_size_{size}.png")

