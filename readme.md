<!-- PARTICLE TEMPERATURE -->
Raihan Alfaridzi - 31 August 2022

Python function for calculating local temperature of an atom, using ovito library. The local temeperature of an atom is defined as the average kinetic energy inside a sphere of radius $R_{cut}$ as

$$T_a=\frac{1}{3 N k_B}* \sum_{i=1}^{N} m_i v_i^2 = 1 $$

where $N$ is the number of neighboring atoms inside the sphere, k_B is Boltzmann's constant, and m_i is the mass of the corresponding atom. The thermal velocit, $v_i$, of atom $i$ is defined as the atom velocity subtracted by the center-of-mass velocity of the sphere.
