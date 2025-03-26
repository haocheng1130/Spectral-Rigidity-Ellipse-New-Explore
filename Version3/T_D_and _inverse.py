from mpmath import mp, mpf, sin, pi, nstr
from sympy.functions.combinatorial.numbers import mobius
import pandas as pd
import pickle

# Set desired precision (e.g., 50 decimal places)
mp.dps = 100

# WARNING: Using a matrix size of 10000x10000 with arbitrary-precision arithmetic
# via pure Python loops will be extremely slow. For testing, use a smaller size.
size = 300

# ------------------------------------------------------------------------------------------------------------
# Step 1: Implement the matrix M and its inverse.
# ------------------------------------------------------------------------------------------------------------
# M is defined by: M[q-1][j-1] = 1 if j is divisible by q, else 0.

# Initialize M as a list of lists (integer values)
M = [[0 for _ in range(size)] for _ in range(size)]
for q in range(1, size + 1):  # 1-based indexing for q
    for j in range(1, size + 1):  # 1-based indexing for j
        if j % q == 0:
            M[q - 1][j - 1] = 1

# M_inverse is defined using the Möbius function:
# If q divides j, then M_inverse[q-1][j-1] = mobius(j // q), else 0.
M_inverse = [[0 for _ in range(size)] for _ in range(size)]
for q in range(1, size + 1):
    for j in range(1, size + 1):
        if j % q == 0:
            M_inverse[q - 1][j - 1] = mobius(j // q)

# ------------------------------------------------------------------------------------------------------------
# Step 2: Implement the matrix D and its inverse.
# ------------------------------------------------------------------------------------------------------------
# Define the sinc function using mpmath for high-precision arithmetic.
def sinc(x):
    # Return sin(x)/x with high precision; define sinc(0)=1.
    return sin(x) / x if x != 0 else mp.mpf(1)

# Create D as a diagonal matrix (stored as a list of lists) with high-precision numbers.
D = [[mp.mpf(0) for _ in range(size)] for _ in range(size)]
D[0][0] = mp.mpf(1)
for q in range(2, size + 1):
    value = pi / mpf(q)  # Convert integer q to high-precision mpf
    D[q - 1][q - 1] = sinc(value)

# Create D_inverse as a diagonal matrix.
D_inverse = [[mp.mpf(0) for _ in range(size)] for _ in range(size)]
D_inverse[0][0] = mp.mpf(1)
for q in range(2, size + 1):
    value = pi / mpf(q)
    D_inverse[q - 1][q - 1] = 1 / sinc(value)

# ------------------------------------------------------------------------------------------------------------
# Step 3: Implement matrix multiplication and compute T_D and T_D_inverse.
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

# Compute T_D = D * M
T_D = mat_mult(D, M)
# Compute T_D_inverse = M_inverse * D_inverse
T_D_inverse = mat_mult(M_inverse, D_inverse)

# ------------------------------------------------------------------------------------------------------------
# Main execution block: Save matrices to CSV files.
# ------------------------------------------------------------------------------------------------------------
def save_matrix_to_csv(matrix_data, filename):
    """
    Save a matrix (list of lists) to a CSV file.
    Each element is converted to a string with a specified number of significant digits.
    """
    with open(filename, "w") as f:
        for row in matrix_data:
            # Convert each element to a string with 10 significant digits
            row_str = ",".join(nstr(x, 10) for x in row)
            f.write(row_str + "\n")

if __name__ == "__main__":

    with open("matrix_T_D.pkl", "wb") as f:
        pickle.dump(T_D, f)
    print("T_D matrix saved to matrix_T_D.pkl")
    
    with open("matrix_T_D_inverse.pkl", "wb") as f:
        pickle.dump(T_D_inverse, f)
    print("T_D_inverse matrix saved to matrix_T_D_inverse.pkl")