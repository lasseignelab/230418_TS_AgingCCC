#!/bin/bash

#SBATCH --job-name=clustering
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=65000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

# set variables
wd="/data/user/tsoelter/projects/230418_TS_AgingCCC"

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER=$USER

# file path to functions script
INPUT="${wd}/src/functions_CCCin3xTgAD.R"

# file path to integrated seurat object
INPUT2="${wd}/data/seurat/integrated_seurat.rds"

# file path to project directory to save outputs and objects
INPUT3="${wd}"

# execute docker
singularity exec \
--cleanenv \
--containall \
-B ${wd} \
${wd}/bin/docker/rstudio_aging_ccc_1.0.0.sif \
Rscript --vanilla ${wd}/src/seurat_preprocessing/02_3xtgad_clustering.R ${INPUT} ${INPUT2} ${INPUT3}
