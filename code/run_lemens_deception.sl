#!/bin/bash

#SBATCH --job-name=lemens
#SBATCH --output=lemens-%A.%a.out # stdout file
#SBATCH --error=lemens-%A.%a.err  # stderr file
#SBATCH --nodes=1
#SBATCH -p short
#SBATCH --time=01:00:00
#SBATCH -A ihesd


# ------------------------------------------------------------------------------
# README:
#
# This script runs the LEMENS workflow to generate drought statistics from
# xanthos using fldgen climate fields that have been downscaled with an2month.
#
#  If you want to run 10,000 realizations of GCM-RCP, and want to submit
#  a job running 100 realizations per node using 100 nodes via a SLURM_ARRAY_TASK_ID
#  you would set 'n_realizations' to 10000 and 'batch_size' to 100 and use
#  '--array=0-99' to spread over 100 nodes.
#
# Example:
#
# sbatch --array=0-99 run_lemens_deception.sl GFDL-ESM2M rcp26 10000 100
#
# ------------------------------------------------------------------------------


# load modules
module purge
source /etc/profile.d/modules.sh >& /dev/null
module load gcc/8.1.0
module load python/miniconda3.9
source /share/apps/python/miniconda3.9/etc/profile.d/conda.sh
module load R/3.4.0

python -m rpy2.situation

START=`date +%s`

python /people/d3y010/projects/lemens/code/lemens_deception.py $SLURM_ARRAY_TASK_ID $1 $2 $3 $4

END=`date +%s`

RUNTIME=$(($END-$START))

echo $RUNTIME
