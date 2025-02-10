from scipy import special
from scipy import integrate
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import time
import timeit
import pickle

# ------------------------------------------------------------------------------------------------------------
# Step1: Follow the paper and Shanza's code,
# we write functions to compute each entry for matrix T with a given eccentricity.
# ------------------------------------------------------------------------------------------------------------

maxq = 10000 # initialize max period number

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

def sinphi_lst(q, str_e): # one possibility to speed this up is to derive an analytic formula for sinphi since we are dealing with ellipses!
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
    sinphi_list_q = sinphi_lst(q,str_e) # this should be cached!!!!!!

    ## This is quite inefficient; one can do this step using numpy probably, or at least generators.
    sum = 0
    for k in range(q):
        sinphi = sinphi_list_q[k]

        laz_k = lazutkin_coordinate_analytic(col_amp[k],e)

        mu_k = mu_analytic(col_amp[k],e)

        sum_exp = sinphi*(np.cos(2*np.pi*j*laz_k)/mu_k)
        sum += sum_exp
    return sum

# ------------------------------------------------------------------------------------------------------------
# Step2: We write functions to compute the reduced matrix T for a given eccentricity.
# ------------------------------------------------------------------------------------------------------------

def lambda_marvizi_melrose(j,str_e,e):
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
    accord=0.01 ## maybe we can use a smaller number !!!!!!!!!!
    start=time.time()
    while continuing:
        q=q+1;
        mmc_new=(q**2 * T_of_q_j(q,j,str_e,e));
        continuing=(np.abs((mmc-mmc_new)/mmc) > accord) # caution: this may throw a div/0 !
        mmc = mmc_new
    elapsed = time.time() - start
    q_trivial = 999
    start=time.time()
    mmc_trivial = (q_trivial**2 * T_of_q_j(q_trivial,j,str_e,e));
    elapsed_trivial = time.time() - start
    print ("Marvizi Melrose coefficient:",j," ",q,"   -   ",mmc_new, mmc_trivial, elapsed, elapsed_trivial)
    return mmc_new
# We probably need to justify the choice of magic_q


# We need to understand the accuracy of the above method; try and plot the values of mmc as q increases
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
# We probably need to justify the choice of magic_q

def reduced_T_qj_matrix(max_q, max_j, str_e, e, lambda_MM):
    """
    Computes the reduced matrix tilde{T}_{q,j} for given parameters.

    Parameters:
        max_q: int
            Maximum number of rows (q) to compute.
        max_j: int
            Maximum number of columns (j) to compute.
        str_e: str
            String representation of the eccentricity (used in semi-axes).
        e: float
            Eccentricity.
        lambda_MM: list
            Precomputed Marvizi-Melrose coefficients for each j.

    Returns:
        reduced_matrix: np.ndarray
            The reduced matrix \tilde{T}_{q,j} of size (max_q, max_j).
    """
    # Initialize the reduced matrix
    reduced_matrix = np.zeros((max_q, max_j))

    for q in range(1, max_q + 1):  # Loop over rows (period q)
        print(f"row {q}")
        for j in range(1, max_j + 1):  # Loop over columns (j)
            # Compute T_{q,j} using the provided function
            T_qj = T_of_q_j(q, j, str_e, e)

            # Fetch \kappa_j = lambda_MM[j-1] (ensure j-1 is valid)
            kappa_j = lambda_MM[j - 1] if j - 1 < len(lambda_MM) else 0

            # Compute \tilde{T}_{q,j}
            reduced_matrix[q - 1, j - 1] = T_qj - kappa_j / (q**2)

    return reduced_matrix



# ------------------------------------------------------------------------------------------------------------
# Step3: Now for a list of eccentricities, I will compute its associated reduceed_T_q_j.
# I will try to store the results in pickle.
# ------------------------------------------------------------------------------------------------------------


########## Some considerations in the interest of computational time ###########################

magic_j = 750  # Arbitrary choice for truncating coefficients
gamma = 3.5  # For faster decay
arbitrary_accuracy = 100  # Controls accuracy of the computation

sampled_e = [0.25]  # Eccentricities to consider,
lambda_MM_dict = {}  # Store Marvizi-Melrose coefficients for each eccentricity
reduced_matrices = {}  # Store reduced matrices for each eccentricity

max_q, max_j = 300, 300 # We can vary these numberes


for e in sampled_e:
    str_e = '%.2f' % e
    collision_pts = pd.read_csv(f"all_periods_{str_e}e_col_amplitudes.txt", sep='\t')
    collision_pts.drop("Unnamed: 0", axis=1, inplace=True)
    print("Collision points loaded")
    test_lambda_marvizi_melrose(10,str_e,e)

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
        if (j % 2 == 0) and (abs(l) < 1e-17): #stopping the computation after the coefficients for even j are close enough to 0
            break
        lambda_MM.append(l)
    print(f"Cached Marvizi-Melrose coefficients for eccentricity {e}; number of MM coefficients cached:{len(lambda_MM)}")
    for j in np.arange(len(lambda_MM) + 1, maxq * arbitrary_accuracy): #at this point the terms were close enough to 0, hence no need to compute
        lambda_MM.append(0)

    lambda_MM_dict[e] = lambda_MM  # Store coefficients

# Second loop: Compute reduced matrices
for e in sampled_e:
    str_e = '%.2f' % e
    print(f"Processing reduced matrix for eccentricity {e}")
    reduced_matrix = reduced_T_qj_matrix(max_q, max_j, str_e, e, lambda_MM_dict[e])
    reduced_matrices[e] = reduced_matrix

# Save the reduced matrices to a pickle file
with open("reduced_matrices.pkl", "wb") as file:
    pickle.dump(reduced_matrices, file)
    print("Reduced matrices saved to reduced_matrices.pkl")


# To load the results from the pickle file
with open("reduced_matrices.pkl", "rb") as file:
    loaded_matrices = pickle.load(file)
    print("Reduced matrices loaded from reduced_matrices.pkl")

# Example: Access the loaded matrix for eccentricity 0.1
print(f"\nReduced matrix for eccentricity 0.1 (first 5 rows):")
print(loaded_matrices[0.1][:5, :5])
