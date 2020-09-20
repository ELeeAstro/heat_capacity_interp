# heat_capacity_interp
Mixes and interpolates heat capacities using JANAF tables

A simple code that reads data from JANAF tables given a list of species and their molecular weights.

Then interpolates c_p to a input temperature and returns:
  - R_bar = the mixed specific gas constants
  - cp_bar = the mixed specific heat capacity
  - kappa_prime = R_bar/cp_bar

A sample test program and call structure for Earth air is provided, the code is well commented.

To use, simply compile e.g. 'gfortran cofp_interp_mod.f90' and run 'a.out'.
When coupling to other programs (e.g. GCM or RT models), comment out the example program.

To add a new species, copy the JANAF data to the 'data' folder and remove the 250, 350 and 450 K lines (if the species has them)
for a total of 64 lines in the JANAF file and rename the file 'species_JANAF.txt'.
