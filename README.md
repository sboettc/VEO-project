# VEO-project
Codes and some data for a project on vectorized Extremal Optimization (VEO)

This repository contains demo-codes for the Extremal Optimization (EO) heuristic.

  (1) "demoSK.c" is a monolithic version of tau-EO that demonstrates how I generated the results in Eur. Phys. J. B 46, 501-505 (2005) (also cond-mat/0407130).
  (2) "parallel.c" is based on (1) and demonstrates VEO (ie. it is only a sequential version of VEO, illustrating its power but not its speed).
  
  Both programs have their own, very primitive, RNG and should compile simply with "gcc file.c -o file.x -lm -O2" and run with "./file.x" to start a self-explanatory dialog. (To produce serial data, the MT-RNG in the gnu-gsl library was used.) 
