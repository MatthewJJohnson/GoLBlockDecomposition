#!/bin/sh

# setting the max wallclock time after which the job will be killed; 
# note: the format is hr:min:sec (default set here to 10 mins)
#SBATCH --time=00:10:00


#"np" stands for number of processes.
#this command will run the job on 2 processes. 
mpirun -np 8 ~/GoLBlockDecomposition/GOL 16 4
