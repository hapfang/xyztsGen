# xyztsGen
The geometries supported include:
 1. A single PWR subchannel geometry
 2. Pipe

An additional line is to be inserted before the probe coordinates in xyzts.dat. For example
	7680           2  1.0000000E-06           5          14           1
 7680          : The total number of probes
 2             : The step skip parameter
 1.0000000E-06 : Not used, let it remain 1.0000000E-06
 5             : Not used, let it remain 5
 14            : Total number of variables (14 - single phase; 15 - two-phase option 1; 19 - two-phase option 2)
 1             : The index of current run (used in ensamble simulations)

