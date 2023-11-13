#!/bin/bash
## run the Rscript PANDA.R and schedule this job to SLURM with
## `sbatch 02_PANDA_array.sh`

#SBATCH --job-name=PANDA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jwhitlock@uab.edu #CHANGE THE EMAIL TO YOUR OWN
#
#Number of tasks needed for this job, the partition, and parameters
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=255000
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=largemem #partition info here: https://docs.rc.uab.edu/cheaha/hardware/#partitions
#
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=0-23 #24 cell type inputs total across both AD and WT fofr astros, micros, oligos, OPCs, ex neurons, and in neurons

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

#load modules
module load Singularity/3.5.2-GCC-5.4.0-2.26

#variables
wd="/data/user/jbarham3/230418_TS_AgingCCC/"
src="/data/user/jbarham3/230418_TS_AgingCCC/src/PANDA" #be sure that your subdirectories are structured the same

export SINGULARITYENV_PASSWORD='pass'
export SINGULARITYENV_USER='jbarham3'

#code to execute docker and script for analysis
cd ${wd}

#array file of cell type specific expression inputs for PANDA
SAMPLE_LIST="${wd}/data/PANDA_inputs/PANDA_exp_files_array.txt" #note: make sure path and file name are correct
SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID

# NOTE user must have already pulled from docker and built .sif file with singularity below (jordanwhitlock/setbp1_manuscript_panda_1.0.1)
singularity exec --cleanenv --containall -B ${wd} ${wd}/bin/docker/setbp1_manuscript_panda_1.0.1_latest.sif Rscript --vanilla ${src}/02_PANDA.R ${INPUT} # here vanilla ensures only the script is run and environment is kept clean