import numpy as np
import pandas as pd
from sympy.ntheory import mobius


# ------------------------------------------------------
# Step1: We implement the matrix M and its inverse. 
# ------------------------------------------------------


# We want the size to be large, so we set size to be 10000*10000 for now
# we can vary the size later.
size = 10000


# Matrix M
M = np.zeros((size, size), dtype=int)
# Populate the matrix based on the divisibility condition
for q in range(1, size + 1):  # 1-based indexing for q
    for j in range(1, size + 1):  # 1-based indexing for j
        if j % q == 0:  # Check if j is divisible by q
            M[q - 1, j - 1] = 1  # Set M[q-1, j-1] = 1 (convert to 0-based indexing)

# Matrix M_inverse
M_inverse = np.zeros((size, size), dtype=int) # Initialize the matrix with zeros

# Populate the matrix based on the divisibility condition and Möbius function
for q in range(1, size + 1):  # Loop over q
    for j in range(1, size + 1):  # Loop over j
        if j % q == 0:  # Check if q divides j
            M_inverse[q - 1, j - 1] = mobius(j // q)  # Compute Möbius function for q/j



# ------------------------------------------------------
# Step2: Implement the matrix D and its inverse. 
# ------------------------------------------------------
            

# Define the sinc function
def sinc(x):
    return np.sin(x) / x if x != 0 else 1.0  # Handle x=0 case to return 1.0

# Matrix D
D = np.zeros((size, size), dtype=float)
# Populate the diagonal entries, notive we need to set the first entry to be 1.
D[0, 0] = 1
for q in range(2, size + 1):
    value = np.pi / q
    D[q - 1, q - 1] = sinc(value)


# Matrix D_inverse
D_inverse = np.zeros((size, size), dtype=float)
# Populate the diagonal entries
D_inverse[0, 0] = 1
for q in range(2, size + 1):
    value = np.pi / q
    D_inverse[q - 1, q - 1] = 1 / sinc(value)



# ------------------------------------------------------
# Step3: Implement the matrix T_D and its inverse.
# ------------------------------------------------------

# Compute T_D = D * M
T_D = np.dot(D, M)
# Compute T_D_inverse = D * M_inverse
T_D_inverse = np.dot(M_inverse, D_inverse) 





# ------------------------------------------------------
# Main execution block for yield curve analysis.
# ------------------------------------------------------
if __name__ == "__main__":
    # Save T_D to a file
    output_file = "matrix_T_D_10000x10000.csv"
    np.savetxt(output_file, T_D, fmt='%.10f', delimiter=',')
    print(f"T_D matrix saved to {output_file}")
    # Save T_D_inverse to a file
    output_file = "matrix_T_D_inverse_10000x10000.csv"
    np.savetxt(output_file, T_D_inverse, fmt='%.10f', delimiter=',')
    print(f"T_D_inverse matrix saved to {output_file}")

# uncomment below for test:
# # We check the first 5 rows and columns of matrix D
# print(D[:5, :5])
# # We check the first 5 rows and columns of matrix D_inverse
# print(D_inverse[:5, :5])