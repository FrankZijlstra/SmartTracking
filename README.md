# SMART Tracking
Reconstruction code for SMART device tracking for MR-guided interventions.
This code accompanies the paper 'SMART tracking: Simultaneous anatomical imaging and real-time passive device tracking for MR-guided interventions'. Please see the paper for details on the method.

## Running
Run download_data.m followed by run_tracking.m in Matlab. At the top of run_tracking
you can modify the dataset name and number of objects to be tracked (1 for needle
datasets, 5 for spheres datasets).

## Acknowledgements
This code uses the NUFFT library by Jeffrey Fessler, part of the Michigan Image
Reconstruction Toolbox, available at https://web.eecs.umich.edu/~fessler/code/index.html.
The library is included in the utils/nufft directory.

This code uses the Munkres algorithm as implemented by Yi Cao, avaiable at
https://nl.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm.
The function is included in the utils directory.
