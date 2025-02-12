from scipy import special, integrate
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import time
import timeit
import pickle
from functools import lru_cache

# ------------------------------------------------------------------------------------------------------------
# Step0: Load the ellipse semi-axis data
# ------------------------------------------------------------------------------------------------------------
semi_axes = pd.read_csv("e_and_semi_axes.txt", sep='\t')
semi_axes.drop("Unnamed: 0", axis=1, inplace=True)

# ------------------------------------------------------------------------------------------------------------
# Step1: Define functions to compute the entries of the T matrix with caching to avoid redundant computations
# ------------------------------------------------------------------------------------------------------------

# Note: Since collision_amplitude and collision_period depend on the global variable collision_pts 
# (which is reloaded for each eccentricity), we include the eccentricity string (str_e) as a parameter 
# so that results for the same eccentricity can be cached.
@lru_cache(maxsize=None)
def collision_amplitude(q, str_e):
    """
    For a given period q, returns a tuple containing the amplitudes (or angles) for q collision points.
    """
    points = collision_pts[str(q).zfill(2)].dropna()
    return tuple(float(points[j]) for j in range(q))

@lru_cache(maxsize=None)
def collision_period(q, str_e):
    """
    For a given period q, returns a tuple containing the 2D coordinates of q collision points.
    Using the global semi_axes data and the amplitude, it calculates:
        x = a * sin(amplitude)
        y = -b * cos(amplitude)
    """
    amplitudes = collision_amplitude(q, str_e)
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    return tuple((a * math.sin(amp), -b * math.cos(amp)) for amp in amplitudes)

def find_vector1(u, v):
    """
    Computes the vector from point u to point v and normalizes it.
    """
    v1 = np.subtract(v, u)
    mag_v1 = math.sqrt(v1[0]**2 + v1[1]**2)
    return tuple(v1 / mag_v1)

@lru_cache(maxsize=None)
def find_tangent_vector(x, y, str_e):
    """
    Computes the unit tangent vector of the ellipse at the point (x, y).
    The semi-axis a and b are provided by semi_axes[str_e].
    """
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    # The computed (x_prime, y_prime) is proportional to the tangent vector.
    y_prime = (b * x) / a
    x_prime = -(a * y) / b
    mag = math.sqrt(x_prime**2 + y_prime**2)
    return (x_prime / mag, y_prime / mag)

@lru_cache(maxsize=None)
def sinphi_lst(q, str_e):
    """
    For a given period q, returns a tuple containing sin(phi) for each collision point,
    where phi is the angle between the trajectory direction and the ellipse's tangent.
    """
    pts = collision_period(q, str_e)
    pts_list = list(pts)
    sinphi_values = []
    for idx, pt in enumerate(pts_list):
        # Previous collision point (using periodic indexing)
        prev_pt = pts_list[(idx - 1) % q]
        v1_unit = find_vector1(pt, prev_pt)
        v2_unit = find_tangent_vector(pt[0], pt[1], str_e)
        cosphi = np.dot(v1_unit, v2_unit)
        # Use max to avoid negative values due to numerical errors
        sinphi_values.append(math.sqrt(max(0, 1 - cosphi**2)))
    return tuple(sinphi_values)

@lru_cache(maxsize=None)
def lazutkin_coordinate_analytic(amplitude, e):
    """
    Computes the Lazutkin coordinate (analytic expression) for a given amplitude and eccentricity e.
    """
    return 0.25 * (special.ellipkinc(amplitude, e**2) / special.ellipk(e**2) - 1)

@lru_cache(maxsize=None)
def mu_analytic(amplitude, e):
    """
    Computes the μ factor, which acts as a normalization factor at the collision point.
    """
    return 2 * special.ellipk(e**2) * math.sqrt((1 - e**2) / (1 - e**2 * math.sin(amplitude)**2))

@lru_cache(maxsize=None)
def T_of_q_j(q, j, str_e, e):
    """
    Computes the entry T(q, j) of the T matrix for a given period q and mode j.
    For q = 1, it directly computes using a fixed angle π/2.
    For q > 1, it sums the contributions from each collision point:
        contribution = sin(phi) * (cos(2π * j * Lazutkin_coordinate) / μ)
    """
    if q == 1:
        mu_k = mu_analytic(math.pi/2, e)
        return 1.0 / mu_k
    amps = collision_amplitude(q, str_e)       # Collision point amplitudes (tuple)
    sinphi_list_q = sinphi_lst(q, str_e)         # sin(phi) for each collision point (tuple)
    s = 0.0
    for k in range(q):
        laz_k = lazutkin_coordinate_analytic(amps[k], e)
        mu_k = mu_analytic(amps[k], e)
        s += sinphi_list_q[k] * (math.cos(2 * math.pi * j * laz_k) / mu_k)
    return s

# ------------------------------------------------------------------------------------------------------------
# Step2: Functions to compute Marvizi–Melrose coefficients and construct the reduced T matrix
# ------------------------------------------------------------------------------------------------------------

def lambda_marvizi_melrose(j, str_e, e):
    """
    For a fixed mode j, computes the Marvizi-Melrose coefficient, i.e., the limit of q^2 * T(q,j) as q → ∞.
    Uses a relative error as the convergence criterion.
    """
    if (j % 2 == 1):
        return 0
    q = j
    mmc = q**2 * T_of_q_j(q, j, str_e, e)
    accord = 0.01
    while True:
        q += 1
        mmc_new = q**2 * T_of_q_j(q, j, str_e, e)
        # Use relative error for convergence; note that mmc is assumed nonzero
        if abs((mmc - mmc_new) / mmc) <= accord:
            mmc = mmc_new
            break
        mmc = mmc_new
    print("Marvizi Melrose coefficient:", j, "q =", q, "value =", mmc_new)
    return mmc_new


# Here is a function to test lambda_marvizi_melrose, which is the only problem left.
def test_lambda_marvizi_melrose(j,str_e,e):
    """
    This computes an approximate limit for q→∞ of T_qj
    it computes the magic_q'th element of the vector

    """
    if (j%2==1): # Ellipse has additional symmetry: every odd term is 0
        return 0
    # magic_q=min(j+10,999)
    q=j;

    mmc=(q**2 * T_of_q_j(q,j,str_e,e));
    continuing=True
    accord=0.001 ## maybe we can use a smaller number !!!!!!!!!!
    start=time.time()
    mmc_list=[];
    while continuing:
        q=q+1;
        mmc_new=(q**2 * T_of_q_j(q,j,str_e,e));
        continuing=(np.abs((mmc-mmc_new)/mmc) > accord) # caution: this may throw a div/0 !
        mmc = mmc_new
        mmc_list.append(mmc)

    elapsed = time.time() - start

    q_trivial = 999
    start=time.time()
    mmc_trivial = (q_trivial**2 * T_of_q_j(q_trivial,j,str_e,e));
    elapsed_trivial = time.time() - start
    print(mmc_list)
    plt.plot(mmc_list[-10:])
    plt.show()
    print ("Marvizi Melrose coefficient:",j," ",q,"   -   ",mmc_new, mmc_trivial, elapsed, elapsed_trivial)
    return mmc_new


def reduced_T_qj_matrix(max_q, max_j, str_e, e, lambda_MM):
    """
    Constructs the reduced T matrix, where each entry is:
         T(q,j) - (lambda_MM[j-1] / q^2)
    """
    reduced_matrix = np.zeros((max_q, max_j))
    for q in range(1, max_q + 1):
        print(f"Processing row {q}...")
        for j in range(1, max_j + 1):
            T_qj = T_of_q_j(q, j, str_e, e)
            kappa_j = lambda_MM[j - 1] if j - 1 < len(lambda_MM) else 0
            reduced_matrix[q - 1, j - 1] = T_qj - kappa_j / (q**2)
    return reduced_matrix

# ------------------------------------------------------------------------------------------------------------
# Step3: Main program – Compute the reduced T matrix for a list of eccentricities and store the results using pickle
# ------------------------------------------------------------------------------------------------------------

magic_j = 750      # Upper limit used for computing Marvizi–Melrose coefficients
gamma = 3.5        # (Currently unused, but may be used to accelerate decay)
arbitrary_accuracy = 100  # Controls the computational accuracy

sampled_e = [0.10]  # List of eccentricities to process
lambda_MM_dict = {}   # Dictionary to store Marvizi–Melrose coefficients for each eccentricity
reduced_matrices = {} # Dictionary to store reduced T matrices for each eccentricity

max_q, max_j = 800, 800  # Dimensions of the reduced T matrix

for e in sampled_e:
    str_e = f"{e:.2f}"
    print(f"\n---\nProcessing eccentricity {e}")
    # Load collision points data; this file contains amplitude data for each period.
    collision_pts = pd.read_csv(f"all_periods_{str_e}e_col_amplitudes.txt", sep='\t')
    collision_pts.drop("Unnamed: 0", axis=1, inplace=True)
    print("Collision points loaded")
    
    # Clear caches for functions that use collision_pts data to avoid using outdated cached results
    collision_amplitude.cache_clear()
    collision_period.cache_clear()
    sinphi_lst.cache_clear()
    lazutkin_coordinate_analytic.cache_clear()
    mu_analytic.cache_clear()
    T_of_q_j.cache_clear()
    
    # Compute Marvizi–Melrose coefficients
    lambda_MM = []
    for j in np.arange(1, magic_j):
        l = lambda_marvizi_melrose(j, str_e, e)
        # For even j, if the coefficient is sufficiently small, stop computing further
        if (j % 2 == 0) and (abs(l) < 1e-7): # we can change the number here
            break
        lambda_MM.append(l)
    # Append zeros for modes that are not computed
    for j in range(len(lambda_MM) + 1, maxq * arbitrary_accuracy):
        lambda_MM.append(0)
    
    lambda_MM_dict[e] = lambda_MM
    print(f"Cached Marvizi-Melrose coefficients for eccentricity {e}")
    
    # Compute the reduced T matrix
    reduced_matrix = reduced_T_qj_matrix(max_q, max_j, str_e, e, lambda_MM)
    reduced_matrices[e] = reduced_matrix

# Save the reduced matrices to a pickle file
with open("reduced_matrices.pkl", "wb") as file:
    pickle.dump(reduced_matrices, file)
    print("Reduced matrices saved to reduced_matrices.pkl")

# Load the results from the pickle file for testing
with open("reduced_matrices.pkl", "rb") as file:
    loaded_matrices = pickle.load(file)
    print("Reduced matrices loaded from reduced_matrices.pkl")

# Print the first 5 rows and first 5 columns of the matrix corresponding to eccentricity 0.1
print(f"\nReduced matrix for eccentricity 0.1 (first 5 rows):")
print(loaded_matrices[0.1][:5, :5])
