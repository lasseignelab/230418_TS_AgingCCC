#!/bin/bash
#SBATCH --job-name CellRanger
#SBATCH --error=logs/%j.%N.err.txt
#SBATCH --output=logs/%j.%N.out.txt
#SBATCH --ntasks=2  ## number of tasks (analyses) to run, i.e. # of nodes to request resource from
#SBATCH --cpus-per-task=4  ## the number of threads allocated to each task
#SBATCH --mem-per-cpu=64G   # memory per CPU core
#SBATCH --partition=medium  ## the partition to run in (short == 12h max run time)
#SBATCH --array=0-11 #12 samples


## Load modules
module load Singularity
module load CellRanger/6.1.1

## Variables
WD="/data/user/tsoelter/230418_TS_AgingCCC"

ID_LIST=(`cat ${WD}/src/cellranger/id_list.txt`) 
ID_INPUT=`echo ${ID_LIST[$SLURM_ARRAY_TASK_ID]}` 

SAMPLE_LIST=(`cat ${WD}/src/cellranger/sample_list.txt`) 
SAMPLE_INPUT=`echo ${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}` 

echo $ID_INPUT
echo $SAMPLE_INPUT

mkdir -p ${WD}/data/CellRangerCounts/pre_soupX/

cd ${WD}/data/CellRangerCounts/pre_soupX/

cellranger count --id=${ID_INPUT} \
		 --transcriptome=/data/user/tsoelter/apps/refdata-gex-mm10-2020-A \
                 --include-introns \
                 --fastqs=/data/project/lasseigne_lab/TabeaSoelter/3xTg_snRNAseq/rawData/ \
                 --sample=${SAMPLE_INPUT}
