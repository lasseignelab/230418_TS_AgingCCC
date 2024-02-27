#!/bin/bash
#SBATCH --job-name=PANDA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@uab.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=largemem
#
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-7

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

# load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

# variables
wd="${USER_DATA}/projects/230418_TS_AgingCCC/"

# change working directory
cd ${wd}

# array file of cell type specific expression inputs for PANDA
SAMPLE_LIST="${wd}/data/PANDA_inputs/PANDA_exp_files_array.txt" 
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) 
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}`

# execute array
singularity exec --cleanenv --containall -B ${wd} ${wd}/bin/docker/setbp1_manuscript_panda_1.0.1_latest.sif \
Rscript --vanilla ${wd}/src/gene_targeting/02_PANDA.R \
${INPUT} ${wd}
