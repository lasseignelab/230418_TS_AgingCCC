# enable args
args <- R.utils::commandArgs(trailingOnly = TRUE)

# set seed
set.seed(2178)
print("seed set")

# load in packages
library(netZooR)
library(data.table)

wd <- args[2]

setwd(wd)
getwd()

# load in functions
source(paste0(wd, "src/functions_JW_PANDA.R"))

# load in the input data needed:
motif <- read.table(file = paste0(wd, "data/PANDA_inputs/mm10_TFmotifs.txt"),
                    sep = "\t")
print("motif loaded")

ppi <- read.table(file = paste0(wd, "data/PANDA_inputs/mm10_ppi.txt"),
                  sep = "\t")
print("ppi loaded")

expression <- loadRData(args[1])
print("expression loaded")

#run panda on multi-omic inputs
pandaResults <- makePanda(motif, ppi, expression)
name <- sub(".Rdata", "", basename(args[1]))
save(pandaResults, file = paste0(wd,
                                 "results/intermediate_outputs/07_panda/",
                                 name,
                                 ".Rdata"))
rm(pandaResults)
print(paste0(name, "_PANDA network made and saved."))