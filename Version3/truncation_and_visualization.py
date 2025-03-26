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
# Part 1: Find the matrix we need to consider
# ------------------------------------------------------------------------------------------------------------

with open("matrix_T_D_inverse.pkl", "rb") as f:
    T_D_inverse = pickle.load(f)  

with open("T_matrices.pkl", "rb") as f:
    T_matrices = pickle.load(f)  

T = T_matrices[0.25]



test_matrix = mat_mult(T_D_inverse, T)

# now test_matrix is a list of lists

# ------------------------------------------------------------------------------------------------------------
# Part 2: Plot all the eigenvalues on the complex plane
# ------------------------------------------------------------------------------------------------------------

for size in range(50, 151, 50):
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
    plt.figure()
    plt.scatter(real_parts, imag_parts, c='blue', marker='o', s=10)
    plt.title(f"Eigenvalue Spectrum for {size}x{size} Submatrix")
    plt.xlabel("Real Part")
    plt.ylabel("Imaginary Part")
    plt.grid(True)
    
    # Save the figure to a PNG file
    plt.savefig(f"spectrum_{size}.png")
    plt.close()
    
    print(f"Saved spectrum_{size}.png")



# ------------------------------------------------------------------------------------------------------------
# Part 2: Plot all the norms of eigenvalues
# ------------------------------------------------------------------------------------------------------------

for size in range(50, 151, 50):
    # Extract top-left 'size' x 'size' submatrix from test_matrix
    submatrix = mp.matrix([row[:size] for row in test_matrix[:size]])
    
    # Compute eigenvalues and eigenvectors (only eigenvalues are needed)
    eigenvalues, _ = mp.eig(submatrix)
    
    # Compute the norm (absolute value) of each eigenvalue
    norms = [abs(ev) for ev in eigenvalues]
    
    # Convert the mpmath numbers to Python floats for plotting
    norms = [n for n in norms]
    
    # Plot the norm of each eigenvalue vs. its index
    plt.figure()
    plt.plot(norms, 'o', markersize=3)  # Adjust markersize as needed
    plt.title(f"Eigenvalue Norms for {size}x{size} Submatrix")
    plt.xlabel("Eigenvalue Index")
    plt.ylabel("Eigenvalue Norm")
    plt.grid(True)
    
    # Save the figure to a PNG file
    plt.savefig(f"norms_{size}.png")
    plt.close()
    
    print(f"Saved norms_{size}.png")
