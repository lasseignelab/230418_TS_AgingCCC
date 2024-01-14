#### Clustering of 3xTgAD snRNA-seq data at multiple resolutions
### Author: Tabea M. Soelter
### Date: 2023-11-21

## Goal: Determine most stable resolution before cell type assignment
## Note: This script is only necessary if you are exploring resolutions.
##       - If using the resolution from our manuscript, you may skip this script

## time tracking
ptm <- proc.time()

## set seed
set.seed(42)
print("seed set")

## load packages
library(tidyverse)
library(clustree)
library(Seurat)
print("packages loaded")

## enable usage of args
args <- R.utils::commandArgs(trailingOnly = TRUE)
print("enabled args usage")

## load my functions
source(args[1])
print("functions loaded")

## load integrated and normalized Seurat object
filtered_seurat <- readRDS(args[2])
print("loaded seurat object")

## clustering
# cluster at multiple resolutions
# I will apply clustree to determine the most stable resolution

# set active assay to RNA
DefaultAssay(filtered_seurat) <- "RNA"

# vector of desired resolutions
resolutions <- c(
  0.5, 0.6, 0.7, 0.8, 0.9,
  1.0, 1.1, 1.2, 1.3, 1.4, 1.5
)

# clustering of resolutions
filtered_seurat <- find_clusters(filtered_seurat,
                                 dims = 1:20,
                                 reduction = "harmony",
                                 resolutions = resolutions
)


print("finished clustering")

# Plot resolutions between 0.5 and 0.9
tree <- clustree(filtered_seurat@meta.data, prefix = "RNA_snn_res.0.")

# Resolutions between 1 and 1.5
tree2 <- clustree(filtered_seurat@meta.data, prefix = "RNA_snn_res.1.")

# Create output file path
filepath <- paste0(args[3],
                   "/results/intermediate_outputs/02_seurat/clustree.pdf")

# Plot clustree outputs
pdf(file = filepath)
tree
tree2
dev.off()
print("saved clustree plots")

# save objects
filepath <- paste0(args[3], "/data/seurat/clustered_seurat.rds")
print("saving clustered object")
saveRDS(filtered_seurat, file = filepath)

# session info
sessionInfo()

# R version 4.2.3 (2023-03-15)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so
# 
# locale:
# [1] C
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] SeuratObject_4.1.3 Seurat_4.3.0       clustree_0.5.0     ggraph_2.1.0      
#  [5] lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0      dplyr_1.1.1       
#  [9] purrr_1.0.1        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1      
# [13] ggplot2_3.4.2      tidyverse_2.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16             colorspace_2.1-0       deldir_1.0-6          
#   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3       
#   [7] spatstat.data_3.0-1    leiden_0.4.3           listenv_0.9.0         
#  [10] farver_2.1.1           graphlayouts_0.8.4     ggrepel_0.9.3         
#  [13] fansi_1.0.4            R.methodsS3_1.8.2      codetools_0.2-19      
#  [16] splines_4.2.3          polyclip_1.10-4        jsonlite_1.8.4        
#  [19] ica_1.0-3              cluster_2.1.4          R.oo_1.25.0           
#  [22] png_0.1-8              uwot_0.1.14            ggforce_0.4.1         
#  [25] shiny_1.7.4            sctransform_0.4.1      spatstat.sparse_3.0-1 
#  [28] compiler_4.2.3         httr_1.4.5             backports_1.4.1       
#  [31] Matrix_1.5-4           fastmap_1.1.1          lazyeval_0.2.2        
#  [34] cli_3.6.1              later_1.3.0            tweenr_2.0.2          
#  [37] htmltools_0.5.5        tools_4.2.3            igraph_1.4.2          
#  [40] gtable_0.3.3           glue_1.6.2             RANN_2.6.1            
#  [43] reshape2_1.4.4         rappdirs_0.3.3         Rcpp_1.0.10           
#  [46] scattermore_0.8        vctrs_0.6.2            nlme_3.1-162          
#  [49] spatstat.explore_3.1-0 progressr_0.13.0       lmtest_0.9-40         
#  [52] spatstat.random_3.1-4  globals_0.16.2         timechange_0.2.0      
#  [55] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
#  [58] irlba_2.3.5.1          goftest_1.2-3          future_1.32.0         
#  [61] MASS_7.3-58.3          zoo_1.8-12             scales_1.2.1          
#  [64] tidygraph_1.2.3        hms_1.1.3              promises_1.2.0.1      
#  [67] spatstat.utils_3.0-2   parallel_4.2.3         RColorBrewer_1.1-3    
#  [70] reticulate_1.28        pbapply_1.7-0          gridExtra_2.3         
#  [73] stringi_1.7.12         checkmate_2.1.0        rlang_1.1.0           
#  [76] pkgconfig_2.0.3        matrixStats_0.63.0     lattice_0.21-8        
#  [79] tensor_1.5             ROCR_1.0-11            labeling_0.4.2        
#  [82] patchwork_1.1.2        htmlwidgets_1.6.2      cowplot_1.1.1         
#  [85] tidyselect_1.2.0       here_1.0.1             parallelly_1.35.0     
#  [88] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
#  [91] R6_2.5.1               generics_0.1.3         pillar_1.9.0          
#  [94] withr_2.5.0            fitdistrplus_1.1-8     survival_3.5-5        
#  [97] abind_1.4-5            sp_1.6-0               future.apply_1.10.0   
# [100] KernSmooth_2.23-20     utf8_1.2.3             spatstat.geom_3.1-0   
# [103] plotly_4.10.1          tzdb_0.3.0             viridis_0.6.2         
# [106] grid_4.2.3             data.table_1.14.8      digest_0.6.31         
# [109] xtable_1.8-4           httpuv_1.6.9           R.utils_2.12.2        
# [112] munsell_0.5.0          viridisLite_0.4.1
# 
# time tracking
ptm <- proc.time() - ptm
fptm <- (fptm[3] / 60) / 60
print(paste0("Run time: ", fptm, " hours"))
# 3 hours and 45 minutes

# Reproducibility:
# This script was styled and linted.
# Code excluded here, as script is submitted as a bash job.