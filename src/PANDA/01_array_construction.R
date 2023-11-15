# NOTE: be sure to run this script to set up list of files BEFORE constructing networks with '02_PANDA_array.sh' bash script. 
# This script was run inside of JW docker container jordanwhitlock/setbp1_manuscript:1.0.5 (setbp1_manuscript_1.0.5.sif) using Singularity/3.5.2-GCC-5.4.0-2.26 
set.seed(2178)

##Load libraries
library(here)
library(utils)

##Make list of .Rdata expression files needed for array job (generated using `MouseSetbp1_PANDA_expression_05.Rmd`):
files <- list.files(here("data/pseudobulk_split"), full.names = TRUE) # ensure that ALL exppression inputs for PANDA are in the same directory and that nothing else is in the directory.
files #save this as a text file to be the array file

write.table(files, file = here("data/PANDA_inputs/PANDA_exp_files_array.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) # for other data user may want to change file name to reflect data and project
