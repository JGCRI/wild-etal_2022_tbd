#!/bin/bash

#SBATCH --job-name=an2month
#SBATCH --output=an2month-%A.%a.out # stdout file
#SBATCH --error=an2month-%A.%a.err  # stderr file
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p short
#SBATCH --time=03:00:00
#SBATCH -A IHESD
#SBATCH --mail-user=chris.vernon@pnnl.gov
#SBATCH --mail-type=ALL


# ------------------------------------------------------------------------------
# README:
#
# This script runs all combinations of model and scenario to generate an2month
# L3 outputs as an RDS file.  There are 16 total combinations.
#
# Example:
#
# sbatch --array=0-15 run_python.sl
#
# ------------------------------------------------------------------------------


# load modules
module purge
module load gcc/8.1.0
module load gdal/2.3.1
module load python/3.7.2
module load R/3.4.3

# activate virtual environment
source /people/d3y010/virtualenvs/py3.7.2_an2month/bin/activate

START=`date +%s`

python /people/d3y010/an2month/code/python_l3.py $SLURM_ARRAY_TASK_ID

END=`date +%s`

RUNTIME=$(($END-$START))

echo $RUNTIME

