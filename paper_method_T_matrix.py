from mpmath import mp, mpf, sin, cos, sqrt, pi, ellipf, ellipk, nstr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle
from functools import lru_cache

# ------------------------------------------------------------------------------
# Set desired arbitrary precision
# ------------------------------------------------------------------------------
mp.dps = 100

# ------------------------------------------------------------------------------
# Step0: Load the ellipse semi-axis data (computed previously with high precision)
# ------------------------------------------------------------------------------
converters = {f"{i/100:.2f}": mp.mpf for i in range(1, 100)}
semi_axes = pd.read_csv("e_and_semi_axes.txt", sep="\t", index_col=0, converters=converters)

# ------------------------------------------------------------------------------
# Step1: Define functions to compute the entries of the T matrix with caching.
#          (All math functions are replaced by mpmath equivalents for high accuracy.)
# ------------------------------------------------------------------------------

# Global variable 'collision_pts' will be loaded in the main loop.
# We include the eccentricity string 'str_e' as a parameter to cache results per eccentricity.

@lru_cache(maxsize=None)
def collision_amplitude(q, str_e):
    """
    For a given period q, return a tuple of amplitudes (angles) for the q collision points.
    Amplitudes are read from the globally loaded collision_pts DataFrame.
    """
    # Retrieve the column corresponding to period q (with 2-digit zero-filled key)
    points = collision_pts[str(q).zfill(2)].dropna()
    # Convert each amplitude to an mpmath high-precision number
    return tuple(mp.mpf(str(points.iloc[j])) for j in range(q))

@lru_cache(maxsize=None)
def collision_period(q, str_e):
    """
    For a given period q, return a tuple of 2D coordinates for the q collision points.
    Coordinates are computed using:
         x = a * sin(amplitude)
         y = -b * cos(amplitude)
    where a and b are the ellipse semi-axes (loaded from the semi_axes DataFrame).
    """
    amplitudes = collision_amplitude(q, str_e)
    # Convert the stored semi-axis values (strings) to high-precision numbers.
    a = mp.mpf(str(semi_axes[str_e].iloc[0]))
    b = mp.mpf(str(semi_axes[str_e].iloc[1]))
    return tuple((a * sin(amp), -b * cos(amp)) for amp in amplitudes)

def find_vector1(u, v):
    """
    Compute the unit vector pointing from point u to point v.
    """
    diff = (v[0] - u[0], v[1] - u[1])
    mag = sqrt(diff[0]**2 + diff[1]**2)
    return (diff[0] / mag, diff[1] / mag)

@lru_cache(maxsize=None)
def find_tangent_vector(x, y, str_e):
    """
    Compute the unit tangent vector of the ellipse at the point (x, y).
    Uses the ellipse's semi-axes (a, b) from semi_axes[str_e].
    The tangent is computed (up to scaling) by differentiating the ellipse equation.
    """
    a = mp.mpf(str(semi_axes[str_e].iloc[0]))
    b = mp.mpf(str(semi_axes[str_e].iloc[1]))
    # Derivative components (proportional to the tangent vector)
    x_prime = - (a * y) / b
    y_prime = (b * x) / a
    mag = sqrt(x_prime**2 + y_prime**2)
    return (x_prime / mag, y_prime / mag)

@lru_cache(maxsize=None)
def sinphi_lst(q, str_e):
    """
    For a given period q, return a tuple of sin(phi) values at each collision point,
    where phi is the angle between the incoming trajectory and the ellipse's tangent.
    """
    pts = collision_period(q, str_e)
    pts_list = list(pts)
    sinphi_values = []
    for idx, pt in enumerate(pts_list):
        # Use periodic indexing for the previous collision point.
        prev_pt = pts_list[(idx - 1) % q]
        v1_unit = find_vector1(pt, prev_pt)
        v2_unit = find_tangent_vector(pt[0], pt[1], str_e)
        # Compute dot product (cosine of the angle)
        cosphi = v1_unit[0] * v2_unit[0] + v1_unit[1] * v2_unit[1]
        # Ensure nonnegative argument for the square root
        sinphi_values.append(sqrt(max(mpf(0), 1 - cosphi**2)))
    return tuple(sinphi_values)

@lru_cache(maxsize=None)
def lazutkin_coordinate_analytic(amplitude, e):
    """
    Compute the Lazutkin coordinate (analytic expression) for a given amplitude and eccentricity e.
    Uses the incomplete and complete elliptic integrals from mpmath.
    """
    # In mpmath, ellipf(phi, m) is the incomplete elliptic integral F(phi|m).
    F_inc = ellipf(amplitude, e**2)
    K_complete = ellipk(e**2)
    return mpf(1) / 4 * (F_inc / K_complete - 1)

@lru_cache(maxsize=None)
def mu_analytic(amplitude, e):
    """
    Compute the μ factor (normalization factor) at the collision point for a given amplitude and eccentricity e.
    """
    K_complete = ellipk(e**2)
    return 2 * K_complete * sqrt((1 - e**2) / (1 - e**2 * (sin(amplitude)**2)))

@lru_cache(maxsize=None)
def T_of_q_j(q, j, str_e, e):
    """
    Compute the entry T(q, j) of the T matrix for a given period q and mode j.
    For q = 1, a fixed angle π/2 is used.
    For q > 1, contributions from each collision point are summed:
         contribution = sin(phi) * [cos(2π * j * Lazutkin_coordinate) / μ]
    """
    if q == 1:
        mu_k = mu_analytic(pi/2, e)
        return mpf(1) / mu_k
    amps = collision_amplitude(q, str_e)  # tuple of amplitudes
    sinphi_list_q = sinphi_lst(q, str_e)
    s = mpf(0)
    for k in range(q):
        laz_k = lazutkin_coordinate_analytic(amps[k], e)
        mu_k = mu_analytic(amps[k], e)
        s += sinphi_list_q[k] * (cos(2 * pi * j * laz_k) / mu_k)
    return s

# ------------------------------------------------------------------------------
# Step2: Build the T matrix.
# ------------------------------------------------------------------------------

def T_qj_matrix(max_q, max_j, str_e, e):
    """
    Constructs the T matrix. The matrix is stored as a NumPy array of Python objects (mpmath numbers).
    """
    T_matrix = np.empty((max_q, max_j), dtype=object)
    for q in range(1, max_q + 1):
        print(f"Processing row {q}...")
        for j in range(1, max_j + 1):
            T_qj = T_of_q_j(q, j, str_e, e)
            T_matrix[q - 1, j - 1] = T_qj
    return T_matrix

# ------------------------------------------------------------------------------
# Step3: Main program – Compute the T matrix for specified eccentricities
#          and store the results using pickle.
# ------------------------------------------------------------------------------

# 'arbitrary_accuracy' controls how many modes are padded with zeros (if needed)
arbitrary_accuracy = 100  

sampled_e = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]  # List of eccentricities to process
# Add other eccentricities if you want to use them

T_matrices = {} # To store reduced T matrices for each eccentricity

max_q, max_j = 300, 300  # Dimensions of the reduced T matrix

for e in sampled_e:
    str_e = f"{e:.2f}"
    print(f"\n---\nProcessing eccentricity {e}")
    # Load collision points data; this file contains amplitude data for each period.
    converters = {f"{q:02d}": lambda x: mp.mpf(x) if x.strip() != "" else "" for q in range(2, maxq)}
    collision_pts = pd.read_csv(f"all_periods_{str_e}e_col_amplitudes.txt", sep='\t', index_col=0, converters=converters)
    print("Collision points loaded")
    
    # Clear caches for functions that depend on collision_pts to avoid outdated results.
    collision_amplitude.cache_clear()
    collision_period.cache_clear()
    sinphi_lst.cache_clear()
    lazutkin_coordinate_analytic.cache_clear()
    mu_analytic.cache_clear()
    T_of_q_j.cache_clear()
    # Note that all the data are cached for computing the reduced T matrix for each eccentricity
    # Since we decided not to use that matrix, we can ignore these codes

    
    
    # Compute the  T matrix
    T_matrix = T_qj_matrix(max_q, max_j, str_e, e)
    T_matrices[e] = T_matrix

# Save the reduced matrices to a pickle file.
with open("T_matrices.pkl", "wb") as file:
    pickle.dump(T_matrices, file)
    print("matrices saved to T_matrices.pkl")