from scipy import special
from scipy import integrate
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle

# ------------------------------------------------------------------------------------------------------------
# Step1: Same functions to compute each entry for matrix T with a given eccentricity.
# ------------------------------------------------------------------------------------------------------------

maxq = 1000 # initialize max period number

semi_axes = pd.read_csv("e_and_semi_axes.txt", sep = '\t') # the file that contains the semi-axes associated with each eccentricity
semi_axes.drop("Unnamed: 0", axis = 1, inplace = True)

def collision_period(q):
    amplitudes = collision_amplitude(q)
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    result = []
    for j in range(0,q):
        result.append([a*math.sin(amplitudes[j]),-b*math.cos(amplitudes[j])])
    return (result)

def collision_amplitude(q):
    """
    Returns an array of the q collision points relative to period q
    """

    result = []
    points = collision_pts[str(q).zfill(2)].dropna()
    for j in range(0,q):
        result.append(float(points[j]))
    return (result)

################### The following functions find the angle associated with each collision point #################################


def find_vector1(u,v): # u is the earlier point, v is the next point/the point we are looking at
    v1 = np.subtract(v,u)
    mag_v1 = math.sqrt((v1[0]**2)+(v1[1]**2))
    v1_unit = v1 / mag_v1
    return (v1_unit)

def find_tangent_vector(x,y,str_e):
    a = semi_axes[str_e][0]
    b = semi_axes[str_e][1]
    y_prime = (b*x)/a
    x_prime = -(a*y)/b
    mag = np.sqrt(x_prime**2+y_prime**2)
    y_tan = y_prime/mag
    x_tan = x_prime/mag
    return (x_tan, y_tan)

def sinphi_lst(q, str_e):
    '''Returns a list of sin(phi) for all point in a given periodic orbit'''
    sinphi_list = []
    for pt in collision_period(q):
        v1_unit = tuple(find_vector1(pt, collision_period(q)[(collision_period(q).index(pt)-1)%q]))
        v2_unit = find_tangent_vector(pt[0],pt[1],str_e)
        cosphi = np.dot(v1_unit,v2_unit)
        sinphi_list.append(np.sqrt(1-cosphi**2))
    return(sinphi_list)

################## The following functions find the lazutkin coordinate associated with each collision point ###########################

def lazutkin_coordinate_analytic(amplitude,e):
    return 0.25*(special.ellipkinc(amplitude,e**2)/special.ellipk(e**2)-1)

def mu_analytic(amplitude,e):
    return 2*special.ellipk(e**2)*np.sqrt((1-e**2)/(1-e**2*math.sin(amplitude)**2))


############# Each entry of the matrix T ###########################

def T_of_q_j(q,j,str_e,e):
    '''T_q_j is a matrix, by varying q and j, you obtain an entry for that matrx.
    Think of the q as the rows, while the j are the columns. Fixing q and j gives you one entry.
    Fixing q and varying j will give you a row, fixing j and varying q will give you a column.
    Make sure to keep this in mind when you use this function, you will have to make changes to the code accordingly
    before you run this function.'''
    k_sum = []
    if q==1:
        mu_k = mu_analytic(np.pi/2,e)
        return 1./mu_k
    col_pt = collision_period(q)
    col_amp = collision_amplitude(q)
    sinphi_list_q = sinphi_lst(q,str_e)

    for k in range(q):
        sinphi = sinphi_list_q[k]

        laz_k = lazutkin_coordinate_analytic(col_amp[k],e)

        mu_k = mu_analytic(col_amp[k],e)

        sum_exp = sinphi*(np.cos(2*np.pi*j*laz_k)/mu_k)
        k_sum.append(sum_exp)
    return(sum(k_sum))


def lambda_marvizi_melrose(j,str_e,e):
    """
    This computes an approximate limit for q→∞ of T_qj
    it computes the magic_q'th element of the vector

    """
    if (j%2==1): # Ellipse has additional symmetry: every odd term is 0
        return 0
    magic_q=min(j+10,999)
    q=j;
    mmc=(q**2*T_of_q_j(q,j,str_e,e));
    continuing=True
    accord=0.000001
    while continuing:
        q=q+1;
        mmc_new=(q**2*T_of_q_j(q,j,str_e,e));
        continuing=(np.abs(mmc-mmc_new)<accord)
    print ("mmc:",j," ",q," ",magic_q,"   -   ",mmc_new)
    return mmc_new
# We probably need to justify the choice of magic_q





# ------------------------------------------------------------------------------------------------------------
# Step2: Using vectorization to compute the reduced matrix T
# ------------------------------------------------------------------------------------------------------------



# ---------------------------
# Precompute Functions
# ---------------------------

def lazutkin_coordinate_analytic_vec(amplitudes, e):
    """Vectorized version for lazutkin_coordinate_analytic."""
    # ellipkinc and ellipk support array inputs in recent SciPy versions.
    # If not, consider using np.vectorize.
    return 0.25*(special.ellipkinc(amplitudes, e**2)/special.ellipk(e**2)-1)

def mu_analytic_vec(amplitudes, e):
    """Vectorized version for mu_analytic."""
    # Ensure vectorized operations:
    return 2*special.ellipk(e**2)*np.sqrt((1-e**2)/(1-e**2*np.sin(amplitudes)**2))

def find_tangent_vector_vec(x, y, a, b):
    """Vectorized tangent vector computation for arrays x, y."""
    # x and y are arrays
    y_prime = (b*x)/a
    x_prime = -(a*y)/b
    mag = np.sqrt(x_prime**2 + y_prime**2)
    return x_prime/mag, y_prime/mag

def sinphi_lst_vec(collision_points, a, b):
    """
    Vectorized version of sinphi_lst for all q simultaneously.
    collision_points is a list or dict: collision_points[q] = array of shape (q,2).
    a,b are semi-axis lengths.
    """
    # For each q, we want sinphi for q points.
    # sinphi = sqrt(1 - cosphi^2), where cosphi = dot(v1_unit,v2_unit)
    # v1_unit is direction from previous point to current
    # v2_unit is tangent vector at current point

    sinphi_dict = {}

    for q, points in collision_points.items():
        # points shape: (q,2)
        # Compute the vectors from prev to current: v1 = p_(k-1)->p_k
        # np.roll points by one to get previous points:
        prev_points = np.roll(points, 1, axis=0)
        v1 = points - prev_points
        mag_v1 = np.linalg.norm(v1, axis=1, keepdims=True)
        v1_unit = v1 / mag_v1

        # Compute tangent vectors at each point
        x, y = points[:,0], points[:,1]
        x_tan, y_tan = find_tangent_vector_vec(x, y, a, b)
        v2_unit = np.column_stack((x_tan, y_tan))

        # cosphi = dot(v1_unit,v2_unit)
        cosphi = np.sum(v1_unit*v2_unit, axis=1)
        sinphi = np.sqrt(1 - cosphi**2)
        sinphi_dict[q] = sinphi
    return sinphi_dict


# ---------------------------
# Main Vectorization Approach
# ---------------------------

def precompute_for_all_q(max_q, str_e, e, collision_pts, semi_axes):
    """
    Precompute all necessary arrays for T_of_q_j computations:
    - collision amplitudes for each q
    - collision points for each q
    - sinphi for each q
    - lazutkin coordinates and mu for each q
    """

    a, b = semi_axes[str_e][0], semi_axes[str_e][1]

    collision_amplitudes = {}
    collision_points = {}
    sinphi_dict = {}
    laz_arr_dict = {}
    mu_arr_dict = {}

    # Precompute ellipk(e^2) once:
    ellipk_e2 = special.ellipk(e**2)

    for q in range(2, max_q+1):
        q_str = str(q).zfill(2) 
        points_series = collision_pts[q_str].dropna()
        col_amp = points_series.values.astype(float)
    
        # Store the amplitudes
        collision_amplitudes[q] = col_amp
    
        x_vals = a * np.sin(col_amp)
        y_vals = -b * np.cos(col_amp)
        pts = np.column_stack((x_vals, y_vals))
        collision_points[q] = pts


    # Compute laz_arr and mu_arr as before, or do this later if you prefer.


    # Compute sinphi for all q in a vectorized manner
    sinphi_dict = sinphi_lst_vec(collision_points, a, b)

    # Compute laz and mu arrays for all q
    for q in range(2, max_q+1):
        col_amp = collision_amplitudes[q]
        # vectorized computations
        laz_arr = lazutkin_coordinate_analytic_vec(col_amp, e)
        mu_arr = mu_analytic_vec(col_amp, e)
        laz_arr_dict[q] = laz_arr
        mu_arr_dict[q] = mu_arr

    return collision_points, sinphi_dict, laz_arr_dict, mu_arr_dict


def reduced_T_qj_matrix_vectorized(max_q, max_j, str_e, e, lambda_MM, collision_points, sinphi_dict, laz_arr_dict, mu_arr_dict):
    """
    Vectorized version of reduced_T_qj_matrix.
    """
    # Ensure lambda_MM is a numpy array and pad if needed
    lambda_MM = np.asarray(lambda_MM)
    if len(lambda_MM) < max_j:
        # pad with zeros
        lambda_MM = np.pad(lambda_MM, (0, max_j - len(lambda_MM)), mode='constant')

    reduced_matrix = np.zeros((max_q, max_j))
    j_array = np.arange(1, max_j+1)

    # For q=1, T_of_q_j is special:
    mu_k = mu_analytic_vec(np.array([np.pi/2]), e)[0]
    T_of_q1_j = 1./mu_k  # constant for all j
    reduced_matrix[0, :] = T_of_q1_j - lambda_MM/(1**2)

    # For q > 1:
    for q in range(2, max_q+1):
        sinphi = sinphi_dict[q]        # shape: (q,)
        laz_arr = laz_arr_dict[q]      # shape: (q,)
        mu_arr = mu_arr_dict[q]        # shape: (q,)

        # Compute cos terms for all j at once:
        # We want cos(2*pi*j*laz_k) for each j and k
        # laz_arr[:, None] * j_array[None, :] gives a q x max_j array
        argument = 2*np.pi*laz_arr[:, None]*j_array[None, :]
        cos_mat = np.cos(argument)  # q x max_j

        # sinphi, mu_arr must be broadcasted appropriately:
        # T_qj = sum over k of (sinphi_k * cos(...)/mu_k)
        # reshape sinphi and mu_arr for broadcasting:
        numerator = sinphi / mu_arr  # q,
        # Perform elementwise multiplication and sum over k:
        T_qj_vector = np.sum((numerator[:, None] * cos_mat), axis=0)

        reduced_matrix[q-1, :] = T_qj_vector - lambda_MM/(q**2)

    return reduced_matrix



# ------------------------------------------------------------------------------------------------------------
# Step3: Now for a list of eccentricities, I will compute its associated reduceed_T_q_j.
# I will try to store the results in pickle.
# ------------------------------------------------------------------------------------------------------------



# Assuming that you've defined or imported:
# semi_axes (from "e_and_semi_axes.txt")
# lambda_marvizi_melrose(j, str_e, e)
# precompute_for_all_q(max_q, str_e, e, collision_pts, semi_axes)
# reduced_T_qj_matrix_vectorized(max_q, max_j, str_e, e, lambda_MM, collision_points, sinphi_dict, laz_arr_dict, mu_arr_dict)

########## Load Semi-Axes Data ##########
semi_axes = pd.read_csv("e_and_semi_axes.txt", sep='\t')
semi_axes.drop("Unnamed: 0", axis=1, inplace=True)

########## Some considerations in the interest of computational time ###########################

magic_j = 750  # Arbitrary choice for truncating coefficients
gamma = 3.5  # For faster decay
arbitrary_accuracy = 100  # Controls accuracy of the computation

sampled_e = [0.10]  # Eccentricities to consider
lambda_MM_dict = {}  # Store Marvizi-Melrose coefficients for each eccentricity
reduced_matrices = {}  # Store reduced matrices for each eccentricity

max_q, max_j = 9999, 9999

# First loop: Compute Marvizi-Melrose coefficients
for e in sampled_e:
    str_e = '%.2f' % e
    print(f"\n---\nProcessing eccentricity {e}")
    
    # Load collision points
    collision_pts = pd.read_csv(f"all_periods_{str_e}e_col_amplitudes.txt", sep='\t')
    collision_pts.drop("Unnamed: 0", axis=1, inplace=True)
    print("Collision points loaded")
    
    
    # Compute Marvizi-Melrose coefficients
    lambda_MM = []
    for j in np.arange(1, magic_j):
        l = lambda_marvizi_melrose(j, str_e, e)
        if (j % 2 == 0) and (abs(l) < 1e-7):
            break
        lambda_MM.append(l)

    # After computing lambda_MM as in Code1:
    for j in np.arange(len(lambda_MM) + 1, maxq * arbitrary_accuracy):
        lambda_MM.append(0)

    # Now truncate to match max_j since we're vectorizing a (max_q, max_j) matrix:
    lambda_MM = np.array(lambda_MM)[:max_j]


    lambda_MM_dict[e] = lambda_MM  # Store the adjusted coefficients

    print(f"Cached Marvizi-Melrose coefficients for eccentricity {e}")

# Second loop: Compute reduced matrices using vectorization
for e in sampled_e:
    str_e = '%.2f' % e
    print(f"Processing reduced matrix for eccentricity {e}")
    
    # Reload collision_pts for use in precomputation
    collision_pts = pd.read_csv(f"all_periods_{str_e}e_col_amplitudes.txt", sep='\t')
    collision_pts.drop("Unnamed: 0", axis=1, inplace=True)
    
    # Precompute all necessary values for vectorized computations
    collision_points, sinphi_dict, laz_arr_dict, mu_arr_dict = precompute_for_all_q(
        max_q, str_e, e, collision_pts, semi_axes
    )
    
    # Compute the reduced matrix using vectorization
    reduced_matrix = reduced_T_qj_matrix_vectorized(
        max_q, max_j, str_e, e, lambda_MM_dict[e], 
        collision_points, sinphi_dict, laz_arr_dict, mu_arr_dict
    )
    reduced_matrices[e] = reduced_matrix

# Save the reduced matrices to a pickle file
with open("reduced_matrices2.pkl", "wb") as file:
    pickle.dump(reduced_matrices, file)
    print("Reduced matrices saved to reduced_matrices2.pkl")

# To load the results from the pickle file later on
with open("reduced_matrices2.pkl", "rb") as file:
    loaded_matrices = pickle.load(file)
    print("Reduced matrices loaded from reduced_matrices2.pkl")

# Example: Access the loaded matrix for eccentricity 0.1
print(f"\nReduced matrix for eccentricity 0.1 (first 5 rows):")
print(loaded_matrices[0.1][:5, :5])
