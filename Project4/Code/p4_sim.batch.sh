#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --partition=amem
#SBATCH --job-name="project4_sim"
#SBATCH --output="/scratch/alpine/rsummers@xsede.org/project4_sim/Logs/%x_%j.out"
#SBATCH --error="/scratch/alpine/rsummers@xsede.org/project4_sim/Logs/%x_%j.err"
#SBATCH --account=amc-general
#SBATCH --time=10:00:00
#SBATCH --mem=400G
#SBATCH --qos=mem

module load miniforge
mamba activate myenv

project_dir="/scratch/alpine/rsummers@xsede.org/project4_sim"

mkdir -p "${project_dir}/Logs" \
"${project_dir}/DataProcessed" \
"${project_dir}/Figures" \
"${project_dir}/Tables"

# directory where run_all script is
cd "${project_dir}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

which Rscript
Rscript -e 'print(.libPaths())'

Rscript run_all.R
