# PARTICLE LOCAL TEMPERATURE FINDER 
Raihan Alfaridzi - 31 August 2022

Python function for calculating local temperature of an atom, using ovito library, for molecular dynamics simulation data. The local temeperature of an atom is defined as the average kinetic energy inside a sphere of radius $R_{cut}$:

$$T_a=\frac{1}{3 N k_B}* \sum_{i=1}^{N} m_i v_i^2 = 1 $$

where $N$ is the number of neighboring atoms inside the sphere, $k_B$ is Boltzmann's constant, and $m_i$ is the mass of the corresponding atom. The thermal velocit, $v_i$, of atom $i$ is defined as the atom velocity subtracted by the center-of-mass velocity of the sphere.


## How to call function

```
local_par_temp(data, saveto, cutoff, units, mass, npy)
```


## Function argument

These are 6 function argument and its explanation :

### data
Directory of input data, that directly pass to the ovito.io import_file function. Reference for importing simulation data can be found in ovito documentation (https://www.ovito.org/docs/current/python/introduction/file_io.html)

### saveto
Directory of where you want to save exported data. There are 2 kind of exported data, which is: xyz data and npy data. npy data is numpy array data which is the array of particel local temperature. xyz data format contains : particle type, X, Y, Z position, temperature. 

### cutoff
Parameter for variable $R_cut$, which is the cutoff radius of neighbor finder. Neighbor of a particle is searched using ovito modifiers, CutoffNeighborFinder.
<br>default = 5.1

### units
Units that are used in simulation. Units type is based on LAMMPS software, that can be found here https://docs.lammps.org/units.html. For example, metal units has mass in grams/mole, and velocity in angstrom/picosecond. This program use only the velocity and mass of particles. All the units are converted into SI units, so the final temperatre result is in Kelvin. 
<br>default = metal

### mass
Input array of particle mass, with the units according to what has been determined previously. Index of array correspond to mass of each atoms respectively. For example if atoms in the simulation are Si; O; H2O; O, then the input array is [28.0855,15.9994,18.0153,15.9994].

### npy
If true, save the file into npy
<br>default = True


## Example
```
data="../10nm/1200/col.*"
saveto="10nm/1200/temp."
Ti=local_par_temp(data, saveto, units="real", mass=[28.0855,15.9994,18.0153,15.9994])
```


