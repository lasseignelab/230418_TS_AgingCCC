README
================
2024-04-04

## Data Directory Structure

As some files in this directory are too large to deposit on GitHub, we
provide a file tree.

Files in this directory are either generated using code from this
project or downloaded from sources specified in our scripts. The
contents of this directory are also deposited on zenodo. Details
(incl. DOIs) can be found in the main repository’s README.

The data directory should include the following files:

    ## .
    ## ├── CellRangerCounts
    ## │   ├── post_soupX
    ## │   │   ├── S01_6m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S02_12m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S03_6m_WT
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S04_12m_WT
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S05_6m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S06_12m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S07_6m_WT
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S08_12m_WT
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S09_6m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S10_12m_AD
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   ├── S11_6m_WT
    ## │   │   │   ├── barcodes.tsv
    ## │   │   │   ├── genes.tsv
    ## │   │   │   └── matrix.mtx
    ## │   │   └── S12_12m_WT
    ## │   │       ├── barcodes.tsv
    ## │   │       ├── genes.tsv
    ## │   │       └── matrix.mtx
    ## │   └── pre_soupX
    ## │       ├── S01_6m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S02_12m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S03_6m_WT
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S04_12m_WT
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S05_6m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S06_12m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S07_6m_WT
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S08_12m_WT
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S09_6m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S10_12m_AD
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       ├── S11_6m_WT
    ## │       │   └── outs
    ## │       │       ├── filtered_feature_bc_matrix.h5
    ## │       │       └── raw_feature_bc_matrix.h5
    ## │       └── S12_12m_WT
    ## │           └── outs
    ## │               ├── filtered_feature_bc_matrix.h5
    ## │               └── raw_feature_bc_matrix.h5
    ## ├── PANDA_inputs
    ## │   ├── PANDA_exp_files_array.txt
    ## │   ├── mm10_TFmotifs.txt
    ## │   └── mm10_ppi.txt
    ## ├── README.Rmd
    ## ├── README.md
    ## ├── ccc
    ## │   ├── 12m_signaling_igraph_objects.rds
    ## │   ├── 6m_signaling_igraph_objects.rds
    ## │   ├── multinichenet_output.rds
    ## │   └── nichenet_v2_prior
    ## │       ├── gr_network_mouse_21122021.rds
    ## │       ├── ligand_target_matrix_nsga2r_final_mouse.rds
    ## │       ├── ligand_tf_matrix_nsga2r_final_mouse.rds
    ## │       ├── lr_network_mouse_21122021.rds
    ## │       ├── signaling_network_mouse_21122021.rds
    ## │       └── weighted_networks_nsga2r_final_mouse.rds
    ## ├── elisa
    ## │   ├── 240319_ELISA_Ab40.csv
    ## │   ├── 240319_ELISA_Ab42.csv
    ## │   └── 240319_ELISA_total_tau.csv
    ## ├── panda
    ## │   ├── excitatory_neurons_AD12.Rdata
    ## │   ├── excitatory_neurons_AD6.Rdata
    ## │   ├── excitatory_neurons_WT12.Rdata
    ## │   ├── excitatory_neurons_WT6.Rdata
    ## │   ├── inhibitory_neurons_AD12.Rdata
    ## │   ├── inhibitory_neurons_AD6.Rdata
    ## │   ├── inhibitory_neurons_WT12.Rdata
    ## │   └── inhibitory_neurons_WT6.Rdata
    ## ├── pseudobulk
    ## │   ├── all_counts_ls.rds
    ## │   ├── astrocytes.rds
    ## │   ├── endothelial_cells.rds
    ## │   ├── ependymal_cells.rds
    ## │   ├── excitatory_neurons.rds
    ## │   ├── fibroblasts.rds
    ## │   ├── inhibitory_neurons.rds
    ## │   ├── meningeal_cells.rds
    ## │   ├── microglia.rds
    ## │   ├── oligodendrocytes.rds
    ## │   ├── opcs.rds
    ## │   ├── pericytes.rds
    ## │   └── rgcs.rds
    ## ├── pseudobulk_split
    ## │   ├── excitatory_neurons_AD12.Rdata
    ## │   ├── excitatory_neurons_AD6.Rdata
    ## │   ├── excitatory_neurons_WT12.Rdata
    ## │   ├── excitatory_neurons_WT6.Rdata
    ## │   ├── inhibitory_neurons_AD12.Rdata
    ## │   ├── inhibitory_neurons_AD6.Rdata
    ## │   ├── inhibitory_neurons_WT12.Rdata
    ## │   └── inhibitory_neurons_WT6.Rdata
    ## └── seurat
    ##     ├── clustered_seurat.rds
    ##     ├── filtered_seurat.rds
    ##     ├── integrated_seurat.rds
    ##     └── processed_seurat.rds
