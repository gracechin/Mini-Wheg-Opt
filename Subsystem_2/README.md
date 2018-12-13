README for the Subsystem_2
==========================
Optimising to minimise the mass of the Mini-Wheg wheel

Scripts
-------

* *latin_hyper_cube.m* - used for sampling
* *area_script.m* - used to check the Rsq normalised model for area
* *mass_script.m* - used to check the Rsq normalised model and real model for mass
* *static_stress_script.m* - used to find and check the Rsq normalised model for static stress
* *drop_stress_script.m* - used to find and check the Rsq normalised model for drop stress  
  
* *std_constraints.m* - normalising constraints
* *optimise.m* - optimises using gradient (with iteration graph) and non-gradient based methods, and outputs results  
  
* *results.csv* - csv of modelling results which is called by other scripts

Execution time
-------
The execution time is approximately 2 seconds

Dependencies
-------
Used in MATLAB R2016b  
Deep Learning Toolbox to use the 'mapstd' function
