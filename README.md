#Parallel Borg training


This repository contains starter code for replication of scaling experiments for the Borg MOEA. To run these, codes, you must also have the master worker and multimaster source code for the Borg MOEA.

These files have been set up and tested to run on The Cube.

To replicate on The Cube:

1. Clone this repository
2. Download and clone the Borg MOEA (http://borgmoea.org) and add to the folder with this repository
3. Load the necessary modules with the following commands:

    module unload gnu8 openmpi3
    module load intel/19.0.2.187

4. Type "make" to compile the code
5. Run the experiment. I have sample DTLZ2 and UF11 SLURM scripts loaded.
5. Calculate hypervolume on runtime files, as demonstrated here: https://waterprogramming.wordpress.com/2019/04/17/performing-random-seed-analysis-and-runtime-diagnostics-with-the-serial-borg-matlab-wrapper/ 


