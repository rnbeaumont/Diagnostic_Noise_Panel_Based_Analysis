#!/bin/bash
#SBATCH --export=ALL    # export environment variables (PBS -V)
#SBATCH -D .    # set working directory to . (PBS -d)
#SBATCH --time=48:0:0   # set walltime (PBS -l walltime)
#SBATCH -p mrcq # queue (PBS -q)
#SBATCH -A Research_Project-MRC158833   # research project to run under
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=count_carriers_Eye
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
##SBATCH --array=1-22	# ARRAY VARIABLE: $SLURM_ARRAY_TASK_ID
##SBATCH --output=%x.%A-%a.out	# FOR ARRAY JOBS
##SBATCH --error=%x.%A-%a.err	# FOR ARRAY JOBS
#SBATCH --mem=50G

PREFIX=Eye
./calculate_revel_quantiles.pl --list ${PREFIX} --out results/${PREFIX}_REVEL_scores
