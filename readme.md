# PARTICLE LOCAL TEMPERATURE FINDER 
Raihan Alfaridzi - 31 August 2022

Python function for calculating local temperature of an atom, using ovito library, for molecular dynamics simulation data. The local temeperature of an atom is defined as the average kinetic energy inside a sphere of radius $R_{cut}$:

$$T_a=\frac{1}{3 N k_B}* \sum_{i=1}^{N} m_i v_i^2 = 1 $$

where $N$ is the number of neighboring atoms inside the sphere, k_B is Boltzmann's constant, and m_i is the mass of the corresponding atom. The thermal velocit, $v_i$, of atom $i$ is defined as the atom velocity subtracted by the center-of-mass velocity of the sphere.


## Function argument

These are 6 function argument and its explanation :

### data
Directory of input data, that directly pass to the ovito.io import_file function. Reference for importing simulation data can be found in ovito documentation (https://www.ovito.org/docs/current/python/introduction/file_io.html)

### saveto
Directory of where you want to save exported data. There are 2 kind of exported data, which is: xyz data and npy data.
