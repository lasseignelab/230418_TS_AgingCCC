#!/bin/bash
#SBATCH --job-name lw_cellranger   ## name that will show up in the queue
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt
#SBATCH --ntasks=1  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=2  ## the number of threads allocated to each task
#SBATCH --mem-per-cpu=12G   # memory per CPU core
#SBATCH --partition=express  ## the partition to run in (short == 12h max run time)
#SBATCH --array=0-10 #11 samples


## Load modules
module load Singularity
module load CellRanger/6.1.1

## Variables
WD="/data/user/lizzyr/code_review/230418_TS_AgingCCC/"
RESULTS="/data/user/lizzyr/code_review/230418_TS_AgingCCC/"

# DATA="/data/project/lasseigne_lab/"
# WD="/data/project/lasseigne_lab/DATASET_dir/George_Lesseigne-selected"
# DOCS="/data/user/tchowton/231027_tc_CitePilot/doc"
# RESULTS="/data/user/tchowton/231027_tc_CitePilot/results/cellranger"

#ID_LIST="${WD}/src/cellranger/id_list.txt"

#ID_INPUT="${ID_LIST[$SLURM_ARRAY_TASK_ID]}"
# ID_INPUT=`echo ${(`cat ${ID_LIST}`)[$SLURM_ARRAY_TASK_ID]}`
ID_LIST=(`cat ${WD}/src/cellranger/id_list.txt`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
ID_INPUT=`echo ${ID_LIST[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID


#ID_INPUT=(`cat ${ID_LIST}`) 
SAMPLE_LIST="${WD}/src/cellranger/sample_list.txt"
#sm_arr=( `cat ${srr_list}|tr "\n" " "`)
# assign sample directory from array and task id
#fq="${sd}/${sm_arr[$SLURM_ARRAY_TASK_ID]}"

# SAMPLE_LIST="${wd}/results/array_inputs/Setbp1_PANDA_files_array.txt" #note: make sure path and file name are correct
# SAMPLE_ARRAY=(`cat ${SAMPLE_LIST}`) # parantheses instruct bash to create a shell array of strings from SAMPLE_LIST
# INPUT=`echo ${SAMPLE_ARRAY[$SLURM_ARRAY_TASK_ID]}` #extracts a single input from the array and prints (using echo) it into INPUT variable, each single input is then assigned an array number by $SLURM_TASK_ID



echo $ID_INPUT
