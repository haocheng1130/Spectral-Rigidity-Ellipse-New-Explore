from mpmath import mp, mpf, sqrt, asin, atan2, ellipe, ellipk, ellipf, ellipfun, pi
import pandas as pd

# Set desired decimal precision (e.g., 100 decimal places)
mp.dps = 100

# ------------------------------------------------------------------------------------------------------------
# Part 1: Create file "e_and_semi_axes.txt"
# For each eccentricity e, compute the semi-major axis a and semi-minor axis b
# for an ellipse of circumference 1.
# ------------------------------------------------------------------------------------------------------------

# Dictionary to hold semi-axis values for each eccentricity.
# Keys are strings (formatted to 2 decimal places).
semi_axis_dict = {}

# Iterate over eccentricities 0.00, 0.01, ..., 0.99
# We can change the range of the eccentricities for testing
for i in range(1, 100):
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
maxq = 500

# Construct a dictionary that maps each eccentricity column name (formatted as "0.00", "0.01", ..., "0.99")
# to the mp.mpf conversion function, ensuring that each value is converted back to a high-precision number.
converters = {f"{i/100:.2f}": mp.mpf for i in range(1, 100)} # be careful about the range of eccentricities!!!

# Read the CSV file "e_and_semi_axes.txt" using tab as the separator.
# Set index_col=0 to treat the first column as the row index (which contains the two rows for a and b).
# Use the converters dictionary to convert each column's string value back into an mp.mpf high-precision number.
semi_axes_df = pd.read_csv("e_and_semi_axes.txt", sep="\t", index_col=0, converters=converters)


# The line of code below is for testing
# semi_axes_df.to_csv("something.txt", sep="\t")

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
    where F is the incomplete elliptic integral of the first kind and
    K is the complete elliptic integral of the first kind.
    """
    a = semi_axes[e][0]
    b = semi_axes[e][1]
    # Compute k_l_sq = (a^2 - b^2) / (a^2 - l^2)
    k_l_sq = (a**2 - b**2) / (a**2 - l**2)
    # Compute the amplitude phi = asin(l/b)
    phi = asin(l / b)
    # Compute the incomplete elliptic integral F(phi, k_l_sq)
    F = ellipf(phi, k_l_sq)
    # Return the rotation number
    return F / (2 * ellipk(k_l_sq))


def find_lambda(w):
    """
    Find an approximate lambda (l) for each eccentricity such that the rotation
    number equals w. A binary search is used for each eccentricity with a tolerance.
    Returns a dictionary mapping eccentricity (string key) to the corresponding lambda.
    """
    l_eccen_dict = {}
    # tol = mp.mpf('1e-10') ### WE MIGHT NEED TO CHANGE THIS!!!
    # If we are working with high precison, this absolute tolerance clearly not doing a good job
    # So instead, we try relative convergence
    tol = mp.mpf('1e-15')   # for typical relative accuracy
    for e in semi_axes:
        a = semi_axes[e][0]
        b = semi_axes[e][1]
        start = mp.mpf(0)
        end = b
        l_val = (start + end) / 2
        while True:
            w_0 = rotation_no(l_val, e)
            if abs(w - w_0) < tol * abs(w):
                l_eccen_dict[e] = l_val
                break
            elif w > w_0:
                start = l_val
                l_val = (start + end) / 2
            else:
                end = l_val
                l_val = (start + end) / 2
            
        # End of binary search for this eccentricity
    return l_eccen_dict

# Precompute lambda values for each period q (as a dictionary).
# period_lambda_dict will map period (as a string) to a dictionary that maps eccentricity to lambda.
period_lambda_dict = {}
for q in range(3, maxq):
    # Compute rotation number w = 1/q (as high-precision number)
    w = mp.mpf(1) / q
    period_lambda_dict[str(q)] = find_lambda(w)
    print(q, "\r")  # simple progress indicator
print("Done finding λ")

def find_collision_pts(e, q):
    """
    For a given eccentricity (e) and period (q), find the collision points.
    The collision points are given by the Jacobi amplitude (phi) computed from the
    Jacobi elliptic functions.
    """
    collisions_dict = {}
    a = semi_axes[e][0]
    b = semi_axes[e][1]
    # Retrieve the corresponding lambda for this eccentricity and period
    l_val = period_lambda_dict[q][e]
    # Compute k_l_sq = (a^2 - b^2) / (a^2 - l^2)
    k_l_sq = (a**2 - b**2) / (a**2 - l_val**2)
    # Compute the complete elliptic integral of the first kind K
    K = ellipk(k_l_sq)
    # For each collision (j=0,...,q-1), compute the collision point.
    for j in range(int(q)):
        # d_l_q divides the interval (4K) evenly among q collisions
        d_l_q = (4 * K) / mp.mpf(q)
        t_j = K + j * d_l_q
        # Compute the Jacobi elliptic functions at t_j with parameter k_l_sq.
        # In mpmath, we use ellipfun: first create the function object for parameter k_l_sq,
        # then evaluate it at t_j. The returned value is a named tuple with attribute 'am' for the amplitude.
        fsn = ellipfun('sn')
        fcn = ellipfun('cn')
        
        # ensure a positive angle, so phi will be in [0, 2pi]
        # I AM NOT SURE HOW TO SET THE RANGE OF phi!!!!!!!!!, because in Shanza's code, phi=am which don't have restriction.
        phi = atan2(fsn(u=t_j, m=k_l_sq),fcn(u=t_j, m=k_l_sq))
        
        if phi<mp.pi/2:
            phi += 2*mp.pi
    
        collisions_dict[str(j).zfill(2)] = phi
    return collisions_dict



# For each eccentricity, compute and store collision points for all periods.
for e in semi_axes:
    # Dictionary to hold collision points for each period for this eccentricity.
    eccen_row_dict = {}
    # Add manually the bouncing ball orbit (period 2) with collision points at pi/2 and 3*pi/2.
    eccen_row_dict["02"] = {"00": mp.pi/2, "01": 3 * mp.pi/2}
    for q in range(3, maxq):
        q_str = str(q)
        eccen_row_dict[str(q).zfill(2)] = find_collision_pts(e, q_str)
        print(q, "\r")  # progress indicator
    # Convert the dictionary to a DataFrame and save to file.
    df_collision = pd.DataFrame(eccen_row_dict)
    filename = f"./all_periods_{e}e_col_amplitudes.txt"
    df_collision.to_csv(filename, sep='\t')
    
print("Done finding collision points")