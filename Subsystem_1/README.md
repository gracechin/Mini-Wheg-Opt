This is the README for the Subsystem_1
=======

Optimisation of the transmission subsystem for minimisation of mass and maximisation of efficiency


Main script 
-------
Final_Optv2.m - Main script to run in this subsystem

Upon completion:

1) xopt : Vector containing all solution variables for fmincon interior point
2) x : Vector containing all solution variables for fmincon global search
3) gax: Vector containing all solution variables for gamultiobj

4) masses: Vector containing optimum values of mass for each of the 3 algorithms
5) effs: Vector containing optimum values of efficiency for each of the 3 algorithms

multiobjfile.m - has multiobjective fitness function for the gamultiobj function call

Execution time
-------
The execution time is approximately 15 seconds

Dependencies
-------
The script requires only MATLAB_R2018A
Deep Learning Toolbox to use the 'mapstd' function
