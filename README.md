# Computational Electromagnetics

This repository contains the series of assignments for the Academic Course "Computational Electromagnetics" taught in the Spring of 2018-2019 in Aristotle University of Thessaloniki - Electrical & Computer Engineering. Both projects are written in Matlab and are focused on **FEM (Finite Element Method)** on **Electrostatic** and **Electromagnetic wave propagation** problems using the **Galerkin** method, weighted residuals formulation.

The below simulations were examined:

<u>Electrostatic:</u> 

- A1 - **Coaxial air-dielectric cable** 
- A2 - **Finite parallel-plate capacitor**

<u>Electromagnetic wave propagation:</u>

- B1 - **Metallic circular waveguide modes (TM/TE)**
- B2 - **Scattering from a Perfect Electromagnetic Conductor cylinder**



## Project A  - Electrostatic Field

The goal of this project is to implement Finite Element Analysis in order to simulate the electric field and the electric potential in a **coaxial air-dielectric cable** (**A1**) as well as in a **finite parallel-plate capacitor (A2)**.  Furthermore other parameters (for example the total energy  stored in the electric field in the given surface) were calculated and compared to the analytical solutions for different number of refinements.



### A1

- Originally , we define the variables that describe the geometry (radius of the interior and external conductor) and construct the geometry description table. We use the **decsg** function to make the decomposed solid geometry matrix, excluding the internal area of the conductor. We make the original grid with **initmesh** function. For 1 refinement we get this geometry:



<p allign = "center">
     <img src="/ProjectA/photos/A1.1.png"width = "70%">
</p>



Using the **table e** and based on the coordinates of the nodes, we locate the nodes corresponding to **Dirichlet conditions** (in this case, the nodes of both conductors - internal and external), thus constructing the **vector node_id**. Also, we construct the **vector X0** that contains the values of the potential at the nodes of the grid, initially only with the known nodes (this vector will be needed for the illustration after the execution of the program). In this vector, nodes of the inner conductor will have a value of V, while those of the outer conductor will have a value of zero. Of course, at the moment the unknown nodes will also have zero potential. 

We create the **index vector** that will contain the number of each unknown nodes (known nodes will have index zero). We initialize the total sparse array and the for loop runs for each element. In each iteration the local tables will be calculated and added to the total one. After that we solve the system SX = B with direct solver (X = S \ B).

We fill the vector X0 (it still has only the values of the known nodes) with unknown values we just found from the solution using the index vector.



- Using the **pdegrad** function, we find and visualize the vector E field.

 <p allign = "center">
     <img src="/ProjectA/photos/A1.2.png"width = "70%">
</p>

In the below table we can compare the results for different number of refinements and write down the Degrees of Freedom in the below table.



| Number of Refinements | Degrees of Freedom |
| :-------------------: | :----------------: |
|           1           |        397         |
|           2           |        1672        |
|           3           |        6861        |



Different solvers of the system vary in time needed. The below table was calculated for **3 refinements** and **1e-6 acccuracy**. 

|             Type of Solver              | Elapsed time (sec) |
| :-------------------------------------: | :----------------: |
| Iterative solver - Biconjugate Gradient |      0.000530      |
|              Direct solver              |      0.037615      |

That is, our system is solved in 71 times faster time with direct relation with iterative solver!



- Based on the result of the field, write a short code to calculate the **total per unit of energy** of the electric field, based on the known relation



<p allign = "center">
     <img src="/ProjectA/photos/energy.png"width = "45%">
</p>


In terms of energy calculation, we know that the field is constant for each triangle and calculated in the center of the element. That means, we have the same number of pairs of Ex and Ey samples as the number of elements. The contribution of each element will be calculated from the square of the electric field at its center on its surface. Once we have found the energy we can solve in terms of its capacity capacitor. The exact value is calculated according to the formula:

```
C = 2πε/ln(b/a)
```

The result of the analytical solution:  Creal=6.6758e-11 F/m

In the below table we can see the relative error of the system, in comparison to the above value of the analytical solution, for different number of refinements.



| Number of Refinements | Relative Error (%) |
| :-------------------: | :----------------: |
|           1           |       2.093        |
|           2           |       1.2939       |
|           3           |       0.6565       |
|           4           |       0.6175       |

We notice a significant improvement with 2 and 3 refinements. Improving the relevant error with 4 refinements is not considered significant.

### A2

Capacitor (infinite length) of parallel plates of finite width w and thickness t has dielectric relative dielectric constant εr between the plates, which are at a distance d, is considered. The calculation area has dimensions Α × Β (Α = Β = 5w) and the capacitor is in the center of it. The potential is V / 2 on the top plate and -V / 2 on the bottom plate. In the external border **Neumann boundary conditions** are considered. We follow the same logic as in A1. The essential difference with respect to A.1 is that for the formation of the vectors node_id and X0, coordinates of the external nodes in table e are  also considered.

- The geometry for 1 refinement: 

<p allign = "center">
     <img src="/ProjectA/photos/A2.1.png"width = "70%">
</p>

- Visualization of potential:

<p allign = "center">
     <img src="/ProjectA/photos/A2.3.png"width = "50%">
</p>


- Visualization of the Electric field:

<p allign = "center">
     <img src="/ProjectA/photos/A2.2.png"width = "70%">
</p>

- Capacity: the only difference with A1 is that we take into consideration only those elements that are located between the 2 plates.

<u>Numerical application</u>: w = 3cm, t = 2mm, d = 1cm, V = 100 Volt, εr = 2.2.

Analytical solution: 
```
C = e_r e_0 w/d
```

C= 5.8436⋅e−11 F/m

In the below table we can see the relative error of the system, in comparison to the above value of the analytical solution, for different number of refinements.

| Number of Refinements | Relative Error (%) |
| :-------------------: | :----------------: |
|           1           |       2.7845       |
|           2           |       1.6578       |
|           3           |       0.7048       |
|           4           |       0.2188       |

We observe very high accuracy here as well. However, the more refinements the slower the code.




## Project B - Electromagnetic Wave Propagation 



### B1

The goal of this project is to illustrate the first 9 **modes (TM or TE) in a metallic circular waveguide**, filled with air. Cutoff frequencies and their relative errors in comparison to the analytical solutions for different number of refinements were also calculated. 

In the case of TE modes, where Hz is unknown, the limit conditions for this component are homogeneous **Neumann** so the nodes of the outer boundary are left as unknown. In contrast, at TE modes the nodes of the outer boundary will have zero values of Ez (**Dirichlet conditions**)



It is an eigenvalue problem, so we need to initialize the two sparse total matrixes of stiffness and mass, S and T, respectively. The local stiffness  matrixes are created completely similar to the case of the electrostatic field. Also, in this case, there is no column vector. 



We solve the generalized problem:


```
(S-k_c^2T)x = 0
```

 using the **eigs** function in order to find the **cutoff wavenumbers** . Use the format [V, D] = eigs (S, T, k, sigma) ;, selecting k = 6 (the first 6 eigenvalues for both TE and TM) and suitable sigma to get the smaller ones eigenvalues.



From the table V containing the eigenvectors we take each of its 6 columns (with a for from 1 to 6) and using it, we fill in the appropriate vector X0 (which will eventually has the values of the unknown field at all nodes) using the vector index.



We visualize the change in the unknown field for each of the first 6 TE modes and first 6 TM modes with **pdeplot**, using 3 refinements.

After running the program a user menu of options for which rhythms will be displayed.



```
If you want to plot the TE modes press 1.
If you want to plot the TM modes press 2.
Make a choice :
```

If the user chooses 1 (TE modes), then, we do not have known nodes therefore the program proceeds directly to the setting of the sparse arrays S, T (After first finding the number of unknown nodes.) In contrast, if the user chooses the 2nd option, the program initializes X0 with zero values for nodes of the external border since it is the TM mode. Afterward, the index table is filled in with the positions of the known nodes of X0 and we can proceed to the core of the algorithm. Once we have found the tables of mass and stiffness we get through eigs the first 6 lower eigenvalues. with the use of 'sm 'in order to get the smaller ones.





- <u>TE01</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TE01.png"width = "50%">
</p>


- <u>TE11</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TE11.png"width = "50%">
</p>


- <u>TE21</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TE21.png"width = "50%">
</p>


- <u>TE31</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TE31.png"width = "50%">
</p>


- <u>TE41</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TE41.png"width = "50%">
</p>


This photo does not represent a mode, so it is ignored.

<p allign = "center">
     <img src="/ProjectB/photos/B1/B1.7.png"width = "50%">
</p>

- <u>TM01</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM01.png"width = "50%">
</p>


- <u>TM02</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM02.png"width = "50%">
</p>


- <u>TM11</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM11.png"width = "50%">
</p>

- <u>TM12</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM12.png"width = "50%">
</p>

- <u>TM21</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM21.png"width = "50%">
</p>

- <u>TM31</u>

<p allign = "center">
     <img src="/ProjectB/photos/B1/TM31.png"width = "50%">
</p>




We can compare our results with the below photo, representing modes of a circular waveguide.

<p allign = "center">
     <img src="/ProjectB/photos/B1/modes.png"width = "70%">
</p>






Now that we have found the cut-off wavenumbers k_c , we can easily find the corresponding cut-off frequencies f_c using the below line:

```
f_c = c*k _c/2π
```


In the below tables we can see the **relative error** of the simulations of the cut-off frequencies (reference the analytical solutions), for different number of refinements.

- <u>TE</u>

| Number of refinements |  TE11  |  TE21  |  TE01  |  TE31  |  TE41  |
| :-------------------: | :----: | :----: | :----: | :----: | :----: |
|           1           | 0.1499 | 0.2746 | 0.4604 | 0.4569 | 0.4326 |
|           2           | 0.0451 | 0.0749 | 0.1099 | 0.1182 | 0.0799 |
|           3           | 0.0188 | 0.0246 | 0.0218 | 0.0330 | 0.0188 |

- <u>TM</u>

| Number of refinements |  TM01  |  ΤΜ11  |  ΤΜ21  | TM02   |  ΤΜ31  |  ΤΜ12  |
| :-------------------: | :----: | :----: | :----: | ------ | :----: | :----: |
|           1           | 0.1501 | 0.3717 | 0.7019 | 0.7646 | 1.0613 | 1.2506 |
|           2           | 0.0323 | 0.0878 | 0.1856 | 0.1934 | 0.2689 | 0.3106 |
|           3           | 0.0026 | 0.0162 | 0.0556 | 0.0495 | 0.0693 | 0.0734 |



The results are excellent both for TE and TM. It seems that we have slightly better accuracy for the first modes. 



### B2

In this project the electromagnetic **scattering from a perfect electromagnetic conductor cylinder**, is considered. The conductor is embedded in the air and a uniform electromagnetic wave on TM mode is approaching it:



<p allign = "center">
     <img src="/ProjectB/photos/B2/E_inc.png"width = "35%">
</p>




A **PML (Perfectly Matched Layer)** was used as an absorbing layer. The field was simulated for different wavelength and different PML thickness. The illustration of the field inside the PML was not included because of the non-existence meaning. The unknown variable of the problem is the z component of the scattered electric field. The distance of the PML (from its four sides) to the cylinder is w ("air") . Below, we can see the geometry and the respective areas.

<p allign = "center">
     <img src="/ProjectB/photos/B2/B2.0.png"width = "50%">
</p>



Completing the geometry we get our original grid with a refinement. (However, we will not need the PML-areas, when we visualize the EM field, for the above reason.)

<p allign = "center">
     <img src="/ProjectB/photos/B2/B2.1.png"width = "50%">
</p>




Since we have defined the field in the known nodes, (those that are located on the boundary of the scatter), we proceed to the core algorithm. Before we build the matrixes of mass and rigidity we need the discretization of the area in order to correctly define the parameters μxx, μyy, εzz in each subarea.

As we create the geometry, Matlab automatically defines the numbering of the regions, which are shown in the image below in blue font.

<p allign = "center">
     <img src="/ProjectB/photos/B2/B2.2.png"width = "40%">
</p>



Parameter a is constant and we choose β so that we have a reflection coefficient of 10^(-6) in the vertical incidence (β depends from the thickness d of the PML).
<p allign = "center">
     <img src="/ProjectB/photos/B2/a_param.png"width = "20%">
</p>
 The formation of the matrixes is very similar to the Electrostatic Problem. Except from the local stiffness matrix, the local mass matrix is calculated in this case too. The total matrix of the system will be a matrix A that will be calculated from the aggregation of the local matrixes. 



<p allign = "center">
     <img src="/ProjectB/photos/B2/A_e.png"width = "20%">
</p>



At the boundary of the scatter, non-homogeneous Dirichlet conditions are applie: 

<p allign = "center">
     <img src="/ProjectB/photos/B2/Dirichlet.png"width = "20%">
</p>

Therefore at each node of the scatter's boundary, the scattering field will be known and equal to the incident. On the contrary, at the external border (termination of the PML) we will apply the most simple boundary conditions, homogeneous Neumann, so the nodes of the external border will be unknown. The arrays node_id, index and X0 will be defined by the scattering boundary nodes.

The script plots different cases for the parameters of diameter and thickness, excluding of course PML-regions. Using 4 refinements mainly for greater accuracy to distinguish differences as thickness of PML changes, we get the below results.  For our arithmetic application a frequency of  300 ΜΗz (λ=1m) and a distance of the conductor from the PML w = λ. Below are given only 2 results for demonstrative purposes, with different d (distance from PML) and radius. More photos with results of the illustrations can be found in the corresponding folder (ProjectB/photos/B2/).



- <u>d = 0.25, a = 0.25</u>

<p allign = "center">
     <img src="/ProjectB/photos/B2/B2.3.png"width = "50%">
</p>


- <u>d = 1, a = 1</u>

<p allign = "center">
     <img src="/ProjectB/photos/B2/B2.7.png"width = "50%">
</p>
