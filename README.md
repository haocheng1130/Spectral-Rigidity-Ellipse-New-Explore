# Spectral-Rigidity-Ellipse-New-Explore
We use the scripts contained herein to provide further numerical evidence for dynamical spectral rigidity among ellipses of various eccentricities. The code contained in these scripts follows the technique implemented in a paper by J. De Simoi, V. Kaloshin, and Q. Wei (a link to the paper is here: [Dynamical_Spectral_Rigidity](https://annals.math.princeton.edu/2017/186-1/p07)) and incorporates portions of code by Shanza Ayub.

There are currently 5 scripts: (as suggested by Professor Jacob De Simoi, I will put everything in py file)

1. `ellipse_axes_and_collisions.py`
2. `T_D_and_inverse.py`
3. `paper_method_reduced_matrix_T.ipynb`
4. `vectorization_reduced_matrix_T.ipynb`
5. `truncation_and_visualization.py`

---

### 1. ellipse_axes_and_collisions.py

**Overview:**  
This script merges the functionality from [Shanza’s two scripts](https://github.com/shanzaayub/Spectral-Rigidity-Ellipse), originally named **e_and_semi_axes_file.py** and **Col_pts_find.py**. The first part calculates the semi-major and semi-minor axes of ellipses (with a fixed circumference of 1) for eccentricities $ e \in [0,1)$, while the second part computes collision points for orbits of different periods $q$ given those same eccentricities.



#### Part 1: Semi-Axes Generation

- **Goal:**  
  Create a file named **e_and_semi_axes.txt** that tabulates the semi-major axis $a$ and semi-minor axis $b$ for each eccentricity $e$ ranging from 0 to 0.99 (with step size 0.01). The circumference is fixed to 1.

- **Key Details:**  
  1. Uses the [Complete Elliptic Integral of the Second Kind](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html) from `scipy.special`.
  2. Calculates:  
     $a = \frac{1}{4 \times E(e^2)}, \quad b = a \sqrt{1 - e^2}$  
     where $E(\cdot)$ is the elliptic integral of the second kind.
  3. Stores each pair $(a,b)$ in a dictionary keyed by the eccentricity $e$ and saves it to **e_and_semi_axes.txt**.



#### Part 2: Collision Points Generation

- **Goal:**  
  For each eccentricity in **e_and_semi_axes.txt**, find the collision (bounce) amplitudes for orbits of period $q$ ranging from 2 up to a user-defined maximum (set here to 10000). Each eccentricity’s collision points are saved in a separate file.

- **Key Steps:**  
  1. **Load Semi-Axes:**  
     Reads **e_and_semi_axes.txt** to retrieve the $(a,b)$ pairs for each eccentricity.
  2. **Define $\lambda$-Finder:**  
     For each period $q$, compute the value $\lambda$ such that the rotation number $\omega = 1/q$ is achieved.
  3. **Collision Points:**  
     Once $\lambda$ is found, use elliptic functions (`scipy.special.ellipj`) to calculate the collision points for each orbit segment (0 through $q-1$).
  4. **Save Results:**  
     For each eccentricity $e$, write all collision points for every period $q$ into **all_periods_\{e\}e_col_amplitudes.txt**.



#### Script Flow Summary

1. **Create e_and_semi_axes.txt:**
   - Loops over $e$ in $\{0.00, 0.01, \ldots, 0.99\}$.
   - Computes the ellipse’s semi-major and semi-minor axes $(a,b)$.
   - Saves the axes in a tab-delimited file **e_and_semi_axes.txt**.

2. **Compute Collision Points for Each Eccentricity and Period:**
   - Reads **e_and_semi_axes.txt**.
   - For each $q$ in $\{2,3,4,\dots,9999\}$:
     - Finds the corresponding $\lambda$ for $\omega = 1/q$.
     - Calculates the collision points using elliptic integrals.
   - Outputs the collision points for each $q$ to **all_periods_\{e\}e_col_amplitudes.txt**.



#### Dependencies

- `numpy`  
- `scipy`  
- `math`  
- `pandas`

Ensure these libraries are installed, for instance via:


#### Usage

1. **Run the script directly** to generate:
   - **e_and_semi_axes.txt** (semi-major, semi-minor axes for each $e$).
   - A series of files **all_periods_\{e\}e_col_amplitudes.txt** for each eccentricity $e$.

2. **Adjust the maximum period** by changing the `maxq` variable (default 10000).  

3. **Check or modify** the step size for eccentricity in the loop (`np.arange(0,1,0.01)`) if you need finer or coarser resolution.


**Note:** This is the **first script** of the repository, merging two of Shanza’s scripts into a single workflow. Make sure to run subsequent scripts in their recommended order if applicable.

---

### 2. T_D_and_inverse.py

**Objective:**  
Implements  six large $\(10000 \times 10000\)$ matrices derived from two core components:  
1. A divisibility-based matrix $M$ (and its inverse using the Möbius function),  
2. A diagonal matrix $D$ (and its inverse using the $\mathrm{sinc}$ function).  

Finally, these matrices are combined to produce $T_D = D \times M$ and $T_D^{-1} = D^{-1} \times M^{-1}$.

**Description:**  
1. **Matrix $M$**: Constructed so that $M_{q,j} = 1$ if and only if $j$ is divisible by $q$.  
2. **Matrix $M^{-1}$**: Uses the Möbius function $\mu\left(\tfrac{j}{q}\right)$ to invert $M$.  
3. **Matrix $D$**: A diagonal matrix whose diagonal entries are defined by $\mathrm{sinc}\left(\tfrac{\pi}{q}\right)$, with the first entry set to 1.  
4. **Matrix $D^{-1}$**: The inverse diagonal matrix whose entries are $1 / \mathrm{sinc}\left(\tfrac{\pi}{q}\right)$.  
5. **Matrix $T_D$**: Defined as $D \times M$.  
6. **Matrix $T_D^{-1}$**: Defined as $D^{-1} \times M^{-1}$.  

**Generates six files:**  
1. **matrix_T_D_10000x10000.csv**  
2. **matrix_T_D_inverse_10000x10000.csv**  


---


### 3. paper_method_reduced_matrix_T.ipynb
 
It strictly follows the method described in the referenced paper (and Shanza’s code) to construct the matrix **T** for a given eccentricity. Each entry of **T** is computed individually, which is why the runtime can become extremely long for larger matrices.

> **Caution**: The code in this script is very slow. It is **not recommended** to run it with large matrix sizes (e.g., $1000 \times 1000$) unless necessary, as the runtime can span many hours.



#### Overview

1. **Step 1**  
   - Implements functions to compute collision points for each billiard trajectory of period $q$.  
   - Defines supporting functions to calculate angles ($\phi$) at each collision, along with Lazutkin coordinates using elliptic integrals.  
   - The main function here is `T_of_q_j()`, which returns a single entry $T_{q,j}$ of the matrix **T**.

2. **Step 2**  
   - Computes the *reduced* matrix $\tilde{T}$ for a given eccentricity.  
   - Incorporates previously calculated Marvizi-Melrose coefficients to adjust each entry $T_{q,j}$ by $\kappa_j / q^2$.  
   - Summarizes everything in the function `reduced_T_qj_matrix()`.

3. **Step 3**  
   - Loops over a list of eccentricities (`sampled_e`) and computes their corresponding $\tilde{T}$ matrices.  
   - Results are saved in a pickle file for later use, so that we do **not** have to re-run these expensive calculations.



#### Key Observations

- A $100 \times 100$ matrix takes around **2 minutes and 30 seconds** to compute.  
- A $1000 \times 1000$ matrix can exceed **900 minutes** (~15 hours).  
- To reduce runtime, the core algorithm should be **optimized** or **parallelized**.



#### Script Details

1. **Collision Point & Lazutkin Coordinates**  
   - `collision_period(q)`: Returns an array of the $q$ collision points $(x, y)$.  
   - `collision_amplitude(q)`: Returns the collision amplitudes needed for various elliptic integral calculations.  
   - `find_vector1(u, v)` / `find_tangent_vector(x, y, str_e)`: Used to determine angles and tangent vectors.  
   - `sinphi_lst(q, str_e)`: Computes $\sin(\phi)$ for each point in a period-q orbit.  
   - `lazutkin_coordinate_analytic()` / `mu_analytic()`: Evaluate elliptic integrals to obtain the Lazutkin coordinate and a factor $\mu$.

2. **Matrix T Construction**  
   - `T_of_q_j(q, j, str_e, e)`:  
     Returns $T_{q,j}$.  
     This summation-based approach depends on the many functions, making it computationally heavy.

3. **Reduced Matrix $\tilde{T}$**  
   - `lambda_marvizi_melrose(j, str_e, e)`:  
     Computes the approximate limit $\kappa_j$ for $T_{q,j}$ as $q \to \infty$.  
   - `reduced_T_qj_matrix(max_q, max_j, str_e, e, lambda_MM)`:  
     Constructs the matrix $\tilde{T}$.

4. **Pickle for Results**  
   - Uses Python’s `pickle` to store computed matrices in a file named `reduced_matrices.pkl`.  
   - You can load these matrices later with:
     ```python
     with open("reduced_matrices.pkl", "rb") as file:
         loaded_matrices = pickle.load(file)
     ```

#### How to Use (only for small size $\tilde{T}$ )

1. **Prepare Data**  
   - Ensure files like `e_and_semi_axes.txt` and `all_periods_XXe_col_amplitudes.txt` (where `XX` is the eccentricity in string format) are in the correct directory.

2. **Adjust Parameters**  
   - `max_q`, `max_j` for matrix size.  
   - `sampled_e` for the list of eccentricities to process.  
   - `magic_j`, `gamma`, `arbitrary_accuracy` for controlling the Marvizi-Melrose calculations.

3. **Run the Script**  
   ```bash
   paper_method_reduced_matrix_T.ipynb


---


### 4. vectorization_reduced_matrix_T.ipynb
 
It follows the same underlying mathematical approach (as laid out in the paper and in Shanza’s code) for computing the matrix **T** of a billiard system with a given eccentricity. **However, the key innovation here is that many of the computations have been vectorized** for efficiency and (potential) speedup over a purely entry-by-entry approach.

> **Note**: While vectorization can significantly reduce runtime, floating-point rounding and the order of arithmetic operations may result in **slightly different numerical values** compared to earlier, purely iterative implementations. This discrepancy typically manifests at the level of $10^{-12}$ to $10^{-15}$.



#### Overview

1. **Step 1**  
   - Implements functions to compute collision points (and associated collision amplitudes) for each period $q$.  
   - Defines helper functions for angles $(\phi$), Lazutkin coordinates, and other geometry-based quantities.

2. **Step 2**  
   - Introduces `lambda_marvizi_melrose(j, str_e, e)`, which provides an approximate limit $\kappa_j$ for $T_{q,j}$ as $q \to \infty$.  
   - The script prints intermediate results for debugging or monitoring convergence.

3. **Step 3**  
   - **Vectorization**: Uses array-based (NumPy) operations to compute $\sin(\phi)$, Lazutkin coordinates, and the entries $T_{q,j}$ in bulk rather than in a nested loop.  
   - Gathers these results into the *reduced* matrix $\tilde{T}$, where  
     \[
       \tilde{T}_{q,j} = T_{q,j} \;-\; \frac{\kappa_j}{q^2}.
     \]

4. **Step 4**  
   - Loops over a list of eccentricities (`sampled_e`) and computes the vectorized $\tilde{T}$ matrices for each.  
   - Saves the results in `reduced_matrices2.pkl`, allowing later retrieval without recomputation.



#### Key Observations

- **Vectorization** significantly restructures the computation, summing many operations into a smaller number of bulk array operations.
- The order and grouping of floating-point arithmetic differ from purely iterative code, so small numerical discrepancies at very high precision are expected.
- With careful implementation, **large** matrix sizes (e.g., $300 \times 300$ or beyond) can be handled more feasibly than before—though performance gains will also depend on hardware and memory constraints.



#### Script Details

1. **Collision Points & Angles**  
   - `collision_amplitude(q)`: Retrieves the $\theta$-like amplitudes for the $q$ collisions from a preloaded `collision_pts` DataFrame.  
   - `collision_period(q)`: Converts these amplitudes to $(x,y)$ coordinates using the ellipse’s semi-axes.  
   - `find_tangent_vector_vec()`: Vectorized calculation of tangent vectors at each collision point.  
   - `sinphi_lst_vec()`: Vectorized computation of $\sin(\phi)$ for an entire array of collision points.

2. **Matrix $T$ Computation**  
   - `T_of_q_j(q, j, str_e, e)`: A reference function (iterative) for computing a single entry $T_{q,j}$. In the vectorized approach, it is replaced by bulk operations:  
     - Summations over \(\sin(\phi) \cdot \cos(2\pi j \cdot \text{lazutkin\_coord}) / \mu\).  
     - Performed across arrays rather than looping over individual indices.

3. **Marvizi–Melrose Coefficients**  
   - `lambda_marvizi_melrose(j, str_e, e)`:  
     Determines $\kappa_j \approx \lim_{q \to \infty} \left(q^2 \cdot T_{q,j}\right)$.  
   - These $\kappa_j$ values are then subtracted off to form $\tilde{T}_{q,j}$.

4. **Reduced Matrix Construction**  
   - `reduced_T_qj_matrix_vectorized()`:  
     Uses the precomputed arrays for $\sin(\phi)$, Lazutkin coordinates, and $\mu$ to fill a 2D NumPy array of size max_q, max_j.  
   - Allows for generating large $\tilde{T}$ matrices more efficiently than a naive nested loop.

5. **Pickle Storage**  
   - Results are saved in `reduced_matrices2.pkl`.  
   - You can reload these arrays later using Python’s `pickle`.



#### Usage and Notes

1. **Input Files**  
   - Requires `e_and_semi_axes.txt`, containing semi-axis lengths for different eccentricities.  
   - Requires collision amplitude files (e.g., `all_periods_XXe_col_amplitudes.txt`, where `XX` is the numeric eccentricity).

2. **Configuration**  
   - `max_q`, `max_j` control the maximum (row, column) size of the matrix.  
   - `sampled_e` is the list of eccentricities to process.  
   - `magic_j`, `arbitrary_accuracy`, etc., are parameters for controlling the Marvizi–Melrose routine.

3. **Run the Script**  
   ```bash
   vectorization_reduced_matrix_T.ipynb

---


### 5. truncation_and_visualization.py

It focuses on **visualizing** the eigenvalues (the spectrum) of the reduced matrices $\tilde{T}$ generated in the previous steps, under various **truncations**. Three distinct visualization approaches are showcased:

1. **Visualization 1** – Plot all eigenvalues for multiple truncations together in a **single** complex-plane plot.  
2. **Visualization 2** – For **each** truncation size, generate **individual** plots.  
3. **Visualization 3** – Create an **animation** over truncation sizes, showing how the eigenvalues in the complex plane evolve as the dimension $N$ increases.


#### Script Breakdown

1. **Loading Data**
   - Reads the precomputed inverse matrix \(`T_D_inverse`\) from a CSV file (e.g., `"matrix_T_D_inverse_10000x10000.csv"`).  
   - Loads the **reduced matrices** $\tilde{T}$ from a pickle file (e.g., `"reduced_matrices.pkl"` or `"reduced_matrices2.pkl"`).  

2. **Truncation Process**
   - A user-defined range of matrix sizes $N$ is specified (e.g., from 10 to 5000 with a step of 10).  
   - For a chosen $N$, the script **truncates** to an $N \times N$ submatrix and computes **its eigenvalues**.  

3. **Visualization 1**
   - **Single Plot** for all chosen truncation sizes.  
   - Each truncation’s eigenvalues are plotted in the **complex plane** (real part vs. imaginary part).  
   - A color map (e.g., `viridis`) distinguishes the different truncation sizes.  
   - Optionally, plots the histogram of eigenvalue norms.

4. **Visualization 2**
   - **One Plot per Truncation Size**.  
   - Iterates through each $N$ in the truncation list, computes the spectrum, and saves a **separate** `.png` file with the eigenvalues in the complex plane.  
   - Incorporates the eccentricity $e$ and the current $N$ in the plot title and filename.

5. **Visualization 3**
   - **Animated** evolution of eigenvalues as $N$ grows.  
   - Uses `matplotlib.animation.FuncAnimation` to generate frames for each truncation size, updating the scatter plot in real time.  
   - Saves the animation as an `.mp4` file (e.g., `"eigenvalues_complex_plane2.mp4"`), allowing a dynamic view of how the spectrum changes with \(N\).

---

#### How the Eigenvalues Are Computed

- First, the script calculates  
  \[
    \text{required\_matrix} = T\_D\_inverse[0{:}N, 0{:}N] \;\times\; \tilde{T}[0{:}N, 0{:}N],
  \]
  effectively combining $\tilde{T}$ with $\mathbf{T}_D^{-1}$ for the chosen truncation dimension $N$.  
- **Eigenvalues** are then found by calling NumPy’s  
  ```python
  eigvals = np.linalg.eigvals(truncated_matrix)



---


## License

[MIT License](LICENSE)

Feel free to open an issue or pull request if you have any improvements or find any issues!
