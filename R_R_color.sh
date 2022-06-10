#!/bin/bash
#SBATCH --qos=regular
#SBATCH --license=SCRATCH
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=R-R-color-try1
#SBATCH --constraint=haswell
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=zg64@duke.edu
#SBATCH --account m1727

export HDF5_USE_FILE_LOCKING=FALSE

source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_13_1/loadLSST.bash
setup lsst_sims
cd $SCRATCH/Rubin-Roman-Redmagic/

python cal_color_parallel.py