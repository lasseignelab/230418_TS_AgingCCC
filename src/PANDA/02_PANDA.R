# NOTE: please create a '.Rprofile' file in your project directory in the same location as your '.Rproj' file. 
# Include a single line with the following code: `R_PROFILE_USER=""`
# This will address any normalizePath warnings looking for the home directory which is not bound within the container
# This script is caled upon and run within the 02_PANDA_array.sh script

#set seed
set.seed(2178)
print("seed set")

#load in package libraries
library(netZooR)
library(data.table)
#library(here)

getwd() #output wd
setwd("/data/user/jbarham3/230418_TS_AgingCCC/")
getwd()
.libPaths() #output libPath

#load in functions
source("/data/user/jbarham3/230418_TS_AgingCCC/src/functions_JW_PANDA.R")

#enable args
args <- R.utils::commandArgs(trailingOnly = TRUE)

#load in the input data needed:
motif <- read.table(file = "/data/user/jbarham3/230418_TS_AgingCCC/data/PANDA_inputs/mm10_TFmotifs.txt", sep = "\t") #load in motif data
print("motif loaded")

ppi <- read.table(file = "/data/user/jbarham3/230418_TS_AgingCCC/data/PANDA_inputs/mm10_ppi.txt", sep = "\t") #load in ppi data
print("ppi loaded")

expression <- loadRData(args[1]) #load expression data from .Rdata in here from $SAMPLE_LIST
print("expression loaded")

#run panda on multi-omic inputs
pandaResults <- makePanda(motif, ppi, expression)
name <- sub(".Rdata", "", basename(args[1]))
save(pandaResults, file = paste0("/data/user/jbarham3/230418_TS_AgingCCC/results/PANDA/", name, ".Rdata"))
rm(pandaResults)
print(paste0(name, "_PANDA network made and saved."))