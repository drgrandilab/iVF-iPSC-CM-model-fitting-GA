#!/bin/bash
#SBATCH --job-name=one_d
#SBATCH --output=res_omp.txt
#SBATCH --error=err_omp.txt
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=26
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# python run_AF_1Hz.py 3


# sh get_all.sh

# matlab -nodisplay -nodesktop -nosplash -r "population_analysis; exit"

matlab -nodisplay -nodesktop -nosplash -r "test_optimisation_ga; exit"
# matlab -nodisplay -nodesktop -nosplash -r "test_optimisation_ga_2; exit"
# python run_ga_main.py
