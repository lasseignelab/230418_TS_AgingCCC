#!/bin/bash
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsoelter@uab.edu
#SBATCH --job-name=CellRanger
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=64000
#SBATCH --nodes=1
#SABTCH --cpus-per-task=4
#SBATCH --time=2-00:00:00
#SBATCH --share
#SBATCH --partition=medium
#SBATCH --error=%A_%a.err.txt
#SBATCH --output=S%A_%a.out.txt
#SBATCH --array=0-11

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################

## Load modules
module load Singularity
module load CellRanger/6.1.1

## Variables
WD="/data/user/tsoelter/projects/230418_TS_AgingCCC"

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
