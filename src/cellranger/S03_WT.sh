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
#SBATCH --error=S03_WT.err.txt
#SBATCH --output=S03_WT.out.txt

########################################
### PUT YOUR COMMANDS BELOW THIS BOX ###
########################################
cd /data/project/lasseigne_lab/TabeaSoelter/3xTg_snRNAseq/CellRanger/outputs/

module load CellRanger/6.1.1

cellranger count --id=S03_6m_WT \
                 --transcriptome=/data/user/tsoelter/apps/refdata-gex-mm10-2020-A \
                 --include-introns \
                 --fastqs=/data/project/lasseigne_lab/TabeaSoelter/3xTg_snRNAseq/rawData/ \
                 --sample=Ta3 \