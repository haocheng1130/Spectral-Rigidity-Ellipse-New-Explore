from mpmath import mp, mpf, sqrt, asin, atan2, ellipe, ellipk, ellipf, ellipfun, pi
import pandas as pd
import concurrent.futures


# Set desired decimal precision (e.g., 50 decimal places)
mp.dps = 50

# ------------------------------------------------------------------------------------------------------------
# Part 1: Create file "e_and_semi_axes.txt"
# For each eccentricity e, compute the semi-major axis a and semi-minor axis b
# for an ellipse of circumference 1.
# ------------------------------------------------------------------------------------------------------------

# Dictionary to hold semi-axis values for each eccentricity.
# Keys are strings (formatted to 2 decimal places).
semi_axis_dict = {}

# Iterate over eccentricities 0.00, 0.01, ..., 0.99
for i in range(0, 100):
    e = mp.mpf(i) / 100  # high-precision eccentricity
    # Compute the semi-major axis "a" using the formula a = 1/(4*E(e^2))
    a = mp.mpf(1) / (4 * ellipe(e**2))
    # Compute the semi-minor axis "b" using b = a * sqrt(1 - e^2)
    b = a * sqrt(1 - e**2)
    key = '{:.2f}'.format(float(e))
    # Store the values as strings (with full precision) so that they are not rounded
    semi_axis_dict[key] = [mp.nstr(a, mp.dps), mp.nstr(b, mp.dps)]

# Create a DataFrame and write the semi-axis data to a tab-separated file.
df_axes = pd.DataFrame(semi_axis_dict)
df_axes.to_csv("e_and_semi_axes.txt", sep='\t')

# ------------------------------------------------------------------------------------------------------------
# Part 2: Create files containing collision points
# for period q in [3, maxq) for eccentricities in [0,1).
# ------------------------------------------------------------------------------------------------------------

# Set maximum period
maxq = 10000

# Read the semi-axis data from file
semi_axes_df = pd.read_csv("e_and_semi_axes.txt", sep='\t')
if "Unnamed: 0" in semi_axes_df.columns:
    semi_axes_df.drop("Unnamed: 0", axis=1, inplace=True)

# Convert the DataFrame into a dictionary where each key (eccentricity)
# maps to a list [a, b] with mpmath high-precision numbers.
semi_axes = {}
for col in semi_axes_df.columns:
    a_str = semi_axes_df[col].iloc[0]
    b_str = semi_axes_df[col].iloc[1]
    a = mp.mpf(a_str)
    b = mp.mpf(b_str)
    semi_axes[col] = [a, b]

# ------------------------------------------------------
# Define functions for collision point computation using mpmath
# ------------------------------------------------------

def rotation_no(l, e):
    """
    For a given lambda (l) and eccentricity (e), compute the rotation number.
    The rotation number is defined as:
        F(asin(l/b), k_l_sq) / (2*K(k_l_sq))
    where k_l_sq = (a^2 - b^2)/(a^2 - l^2)
    """
    a, b = semi_axes[e]
    k_l_sq = (a**2 - b**2) / (a**2 - l**2)
    phi = asin(l / b)
    F_val = mp.ellipf(phi, k_l_sq)
    K_val = mp.ellipk(k_l_sq)
    return F_val / (2 * K_val)

def find_lambda(w):
    """
    Find an approximate lambda (l) for each eccentricity such that the rotation
    number equals w. Uses a binary search with tolerance 1e-10.
    Returns a dictionary mapping eccentricity (string key) to the corresponding lambda.
    """
    l_eccen_dict = {}
    tol = mp.mpf('1e-10')
    for e in semi_axes:
        a, b = semi_axes[e]
        start = mp.mpf(0)
        end = b
        l_val = (start + end) / 2
        while True:
            w_0 = rotation_no(l_val, e)
            if abs(w - w_0) < tol:
                l_eccen_dict[e] = l_val
                break
            elif w > w_0:
                start = l_val
            else:
                end = l_val
            l_val = (start + end) / 2
    return l_eccen_dict

# Precompute lambda values for each period.
# period_lambda_dict maps period (as a string) to a dictionary that maps eccentricity to lambda.
period_lambda_dict = {}
for q in range(3, maxq):
    w = mp.mpf(1) / q
    period_lambda_dict[str(q)] = find_lambda(w)
    print(q, "\r")
print("Done finding λ")

def find_collision_pts(e, q):
    """
    For a given eccentricity (e) and period (q), find the collision points.
    The collision points are given by the Jacobi amplitude (phi) computed via mpmath's ellipfun.
    """
    collisions_dict = {}
    a, b = semi_axes[e]
    l_val = period_lambda_dict[q][e]
    k_l_sq = (a**2 - b**2) / (a**2 - l_val**2)
    K_val = mp.ellipk(k_l_sq)
    for j in range(int(q)):
        d_l_q = (4 * K_val) / mp.mpf(q)
        t_j = K_val + j * d_l_q
        # Create the function objects for sn and cn with parameter k_l_sq.
        fsn = mp.ellipfun('sn')
        fcn = mp.ellipfun('cn')
        # Compute the amplitude using atan2 to mimic a positive angle.
        phi = atan2(fcn(u=t_j, m=k_l_sq), fsn(u=t_j, m=k_l_sq))
        
        collisions_dict[str(j).zfill(2)] = phi
    return collisions_dict

def compute_collision_for_eccentricity(e):
    """
    Compute collision points for all periods for a given eccentricity e.
    Returns a tuple (e, eccen_row_dict) where eccen_row_dict maps period strings to collision dictionaries.
    """
    eccen_row_dict = {}
    # Manually add the bouncing ball orbit (period 2)
    eccen_row_dict["02"] = {"00": mp.pi/2, "01": 3 * mp.pi/2}
    for q in range(3, maxq):
        q_str = str(q)
        eccen_row_dict[str(q).zfill(2)] = find_collision_pts(e, q_str)
        # Optionally, print progress per eccentricity
        print(f"Eccentricity {e}, period {q}")
    return e, eccen_row_dict

# Parallelize over eccentricities using ProcessPoolExecutor.
results = []
with concurrent.futures.ProcessPoolExecutor() as executor:
    # Dispatch computation for each eccentricity.
    futures = {executor.submit(compute_collision_for_eccentricity, e): e for e in semi_axes}
    for future in concurrent.futures.as_completed(futures):
        e, eccen_row_dict = future.result()
        # Convert the dictionary to a DataFrame and save to file.
        df_collision = pd.DataFrame(eccen_row_dict)
        filename = f"./all_periods_{e}e_col_amplitudes.txt"
        df_collision.to_csv(filename, sep='\t')
        print(f"Done for eccentricity {e}")
        
print("Done finding collision points")