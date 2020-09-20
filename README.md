# heat_capacity_interp
Mixes and interpolates heat capacities using JANAF tables

A simple code that reads data from JANAF tables given a list of species and their molecular weights.

Then interpolates c_p to a input temperature and returns 
  - R_bar = the mixed specific gas constants
  - cp_bar = the mixed specific heat capacity
  - kappa_prime = R_bar/cp_bar 
  
A sample test program and call structure for Earth air is provided, the code is well commented.
