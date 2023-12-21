### Functions for CCC in 3xTgAD project 
# Tabea M. Soelter 

## remove_ambientRNA
# A function which removes ambient RNA from h5 files from Cell Ranger (inputs) and generates filtered barcodes, features, and matrix tsv files, which serve as input to generate Seurat objects during preprocessing. 
remove_ambientRNA <- function(inputs, outputs, plots) {
  print("Making list of objects")
  counts_list <- list.dirs(inputs, 
                           full.names = TRUE, 
                           recursive = FALSE)
  pdf(paste0(plots, "rho_density_plots.pdf"))
  for (i in counts_list) {
    set.seed(42)
    # load in cell ranger h5 outputs
    print("Loading cell ranger h5 objects")
    filt_matrix <- Read10X_h5(paste0(i, "/outs/filtered_feature_bc_matrix.h5"))
    raw_matrix <- Read10X_h5(paste0(i, "/outs/raw_feature_bc_matrix.h5"))
    # create seurat object
    print("Making seurat object")
    object <- CreateSeuratObject(counts = filt_matrix)
    # make soup channel object required for soupX
    print("Making soup channel object")
    sco <- SoupChannel(raw_matrix, filt_matrix)
    # get cluster info (soupX requires very basic level clustering info)
    print("Get cluster info")
    object <- SCTransform(object, verbose = FALSE)
    object <- RunPCA(object, approx = FALSE, verbose = FALSE)
    object <- RunUMAP(object, dims = 1:30, verbose = FALSE)
    object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
    object <- FindClusters(object, verbose = FALSE)
    # adding metadata to soup channel object needed for automatic estimation of ambient RNA
    meta <- object@meta.data
    umap <- object@reductions$umap@cell.embeddings
    sco <- setClusters(sco, setNames(meta$seurat_clusters, rownames(meta)))
    sco <- setDR(sco, umap)
    # Analyzing the soup (automatic - function from soupX)
    # Estimates level of ambient RNA in each sample and removes cells with high levels
    print("Profiling the soup")
    sco <- autoEstCont(sco)
    # Create integer matrix
    adjusted_matrix <- adjustCounts(sco, roundToInt = TRUE)
    # Create output directory if needed
    if (!dir.exists(outputs)) dir.create(outputs)
    # Save filtered barcodes, features, and matrix tsv files
    sample_name <- basename(i)
    filename <- paste0(outputs, sample_name)
    print(paste0("Saving filtered objects to: ", filename))
    DropletUtils::write10xCounts(filename, adjusted_matrix) 
  }
  dev.off()
}

## make_seurat_object
# A function which takes a path to sample folders with the three CellRanger output files and creates a merged seurat object
make_seurat_object <- function(path){
  counts_list <- list.dirs(path, 
                           full.names = TRUE, 
                           recursive = FALSE)
  object_list <- vector('list')
  for (i in counts_list){
    counts <- Read10X(i)
    print(i)
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    seurat_object <- CreateSeuratObject(counts = counts, project = sample_name, min.features = 200)
    if (str_sub(sample_name, - 2, - 1) == "AD") {
      seurat_object$orig.ident <- "AD"
      if (str_sub(sample_name, - 4, - 1) == "6mAD") {
        seurat_object$timepoint <- "6m"
      } else {
        seurat_object$timepoint <- "12m"
      }
    } else {
      seurat_object$orig.ident <- "WT"
      if (str_sub(sample_name, - 4, - 3) == "6m") {
        seurat_object$timepoint <- "6m"
      } else {
        seurat_object$timepoint <- "12m"
      }
    }
    object <- seurat_object
    object_list[[i]] <- object
  } 
  print("Making 6m diseased Seurat Object")
  AD_list <- object_list[grepl("6m_AD", names(object_list))]
  for (i in names(AD_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    AD_list[[i]] <- RenameCells(AD_list[[i]],
                                add.cell.id = sample_name)
  }
  diseased_6m <- Merge_Seurat_List(AD_list)
  
  print("Making 12m diseased Seurat Object")
  AD_list <- object_list[grepl("12m_AD", names(object_list))]
  for (i in names(AD_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    AD_list[[i]] <- RenameCells(AD_list[[i]],
                                add.cell.id = sample_name)
  }
  diseased_12m <- Merge_Seurat_List(AD_list)
  
  print("Merging 6m and 12m diseased objects")
  merged_seurat_AD <- merge(x = diseased_6m,
                            y = diseased_12m,
                            add.cell.id = c("6m", "12m"))
  
  print("Making 6m wildtype Seurat Object")
  WT_list <- object_list[grepl("6m_WT", names(object_list))]
  for (i in names(WT_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    WT_list[[i]] <- RenameCells(WT_list[[i]],
                                  add.cell.id = sample_name)
  }
  wildtype_6m <- Merge_Seurat_List(WT_list)
  
  print("Making 12m wildtype Seurat Object")
  WT_list <- object_list[grepl("12m_WT", names(object_list))]
  for (i in names(WT_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    WT_list[[i]] <- RenameCells(WT_list[[i]],
                                  add.cell.id = sample_name)
  }
  wildtype_12m <- Merge_Seurat_List(WT_list)
  
  print("Merging 6m and 12m wildtype objects")
  merged_seurat_WT <- merge(x = wildtype_6m,
                            y = wildtype_12m,
                            add.cell.id = c("6m", "12m"))
  
  print("Making merged Seurat Object")
  merged_seurat <- merge(x = merged_seurat_WT,
                         y = merged_seurat_AD,
                         add.cell.id = c("WT", "AD"))
  return(merged_seurat)
}

## calculate_qc
# A function which calculates quality control metrics for a merged seurat object
calculate_qc <- function(seurat_object){
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, 
                                                  pattern = "^mt-")
  seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
  return(seurat_object)
}

## format_metadata
# A function which extracts and formats metadata of a seurat object
format_metadata <- function(seurat_object){
  metadata <- seurat_object@meta.data
  metadata$cells <- rownames(metadata)
  # Create sample columns -----
  metadata$sample <- NA
  metadata$sample <- ifelse(
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 4, 5) == "6m",
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 7, 9),
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 8, 10)
  )
  
  metadata$orig_sample_id <- ifelse(
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 4, 5) == "6m",
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 7, 13),
    str_sub(merged_seurat@assays[["RNA"]]@data@Dimnames[[2]], 8, 15)
  )
  # Rename columns -----
  metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                         nUMI = nCount_RNA,
                                         nGene = nFeature_RNA)
  return(metadata)
}

## plot_qc
# A function which takes seurat metadata and plots quality control metrics for filtering purposes
plot_qc <- function(metadata) {
  # Visualize the number of cell counts per condition
  number_of_cells <- metadata %>% 
    ggplot(aes(x = seq_folder, fill = seq_folder)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells") 
  # Visualize the number UMIs/transcripts per cell
  number_of_umis <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nUMI, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  # Visualize the distribution of genes detected per cell
  dist_genes_per_cell <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nGene, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  novelty_score <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = seq_folder, fill = seq_folder)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  # Visualize the distribution of mitochondrial gene expression detected per cell
  dist_mito_gex <- metadata %>% 
    ggplot(aes(color = seq_folder, x = mitoRatio, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells   with low numbers of genes/UMIs
  cor <- metadata %>% 
    ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~seq_folder)
  # Plot QC metrics
  plot(number_of_cells) 
  plot(number_of_umis)
  plot(dist_genes_per_cell)
  plot(novelty_score)
  plot(dist_mito_gex)
  plot(cor)
}

## cell_cycle_effects
# A function which calculates and plots the effect of cell cycle on the data using a filtered seurat object as input. It also performs log normalization, scaling, and dimension reduction using PCA
cell_cycle_effects <- function(filtered_seurat, g2m_genes, s_genes){
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m_genes,
                                      s.features = s_genes)
  filtered_seurat <- FindVariableFeatures(filtered_seurat,
                                          selection.method = "vst",
                                          verbose = FALSE)
  # scale data -----
  filtered_seurat <- ScaleData(filtered_seurat)
  # run pca -----
  set.seed(42)
  filtered_seurat <- RunPCA(filtered_seurat, approx = FALSE)
  # plot pca -----
  elbow <- ElbowPlot(filtered_seurat, reduction = "pca", ndims = 50)
  # plot cell cycle scoring -----
  cell_cycle_plot <- DimPlot(filtered_seurat,
                             reduction = "pca",
                             group.by = "Phase",
                             split.by = "Phase")
  plot(cell_cycle_plot)
  plot(elbow)
  return(filtered_seurat)
}

## harmony_integration
# A function which integrates a seurat object using harmony
harmony_integration <- function(seurat_object, dims){
  set.seed(42)
  seurat_object <- RunHarmony(seurat_object,
                              group.by.vars = "sample",
                              reduction = "pca",
                              dims.use = dims, assay.use = "RNA")
  # Here we use all PCs computed from Harmony for UMAP calculation -----
  set.seed(42)
  seurat_object <- RunUMAP(seurat_object, dims = dims, reduction = "harmony", reduction.name = "umap_harmony")
  return(seurat_object)
}

## find_clusters
# A function which finds clusters for a Seurat Object at user determined resolutions
find_clusters <- function(object, dims, reduction, resolutions) {
  # set reduction method to harmony -----
  object <- FindNeighbors(object, 
                          dims = dims, 
                          reduction = reduction)
  # clustering (Leiden aka algorithm 4)
  for (res in resolutions) {
    print(res)
    object <- FindClusters(object,
                           graph.name = "RNA_snn",
                           resolution = res,
                           algorithm = 4,
                           method = "igraph")
  } 
  return(object)
}

## find_markers
# A function allowing for individual FindMarkers analyses for a defined list of clusters
find_markers <- function(object, resolution, identities, value){
  # set resolution ----------
  object <- SetIdent(object, value = resolution)
  # create empty df ----------
  top_markers <- data.frame()
  # iterating through all unidentified clusters ----------
  for (i in identities) {
    # find marker genes ----------
    markers <- FindMarkers(object, ident.1 = i, max.cells.per.ident = 100, logfc.threhold = 0.25, only.pos = TRUE)
    # print out which cluster was completed ----------
    print(paste0("Markers for cluster ", i))
    # filter markers by specified value and adjusted p-value ----------
    markers_filt <- markers %>% top_n(-value, p_val_adj) %>% add_column(cluster = i)
    # bind empty df and filtered markers together ----------
    top_markers <- rbind(top_markers, markers_filt)
  }
  return(top_markers)
}

## make_sce
# A function which requires a Seurat object in order to create single cell experiment object. It also appends the necessary metadata for pseudo-bulking by cell type and sample.
make_sce <- function(object) {
  # raw data
  counts <- object@assays$RNA@counts
  # metadata
  metadata <- object@meta.data
  # add sample_id column as class 'factor'
  metadata$sample_id <- metadata$sample %>% as.factor()
  # add condition information
  metadata$group_id <- metadata$orig.ident
  metadata$group_id <- relevel(metadata$group_id, "WT")
  # add cell type information
  metadata$cluster_id <- factor(object@active.ident)
  # make sce object
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
  
  return(sce)
}

## pseudobulk
# A wrapper function which generates pseudo-bulked data from single cell experiment objects by cell type. The code was adapted form the HBC pseudobulk tutorial.
pseudobulk <- function(sce) {
  # prepare for pseudo-bulking
  cluster_names <- levels(colData(sce)$cluster_id)
  print(paste0(length(cluster_names), " cell types"))
  sample_names <- levels(colData(sce)$sample_id)
  print(paste0(length(sample_names), " samples"))
  groups <- colData(sce)[, c("cluster_id", "sample_id")]
  aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                  groupings = groups, fun = "sum")
  aggr_counts <- t(aggr_counts)
  # Loop over cell types and extract counts (pseudo-bulk)
  counts_ls <- list()
  for (i in 1:length(cluster_names)) {
    column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
    counts_ls[[i]] <- aggr_counts[, column_idx]
    names(counts_ls)[i] <- cluster_names[i]
  }
  return(counts_ls)
}

## cts_metadata
# A wrapper function to create cell-type-specific metadata for the previously pseudo-bulked count data. This code was adapted from the HBC pseudobulk tutorials.
cts_metadata <- function(sce, counts_list) {
  metadata <- colData(sce) %>% 
    as.data.frame() %>% 
    dplyr::select(group_id, sample_id, timepoint)
  metadata <- metadata[!duplicated(metadata), ]
  rownames(metadata) <- metadata$sample_id
  t <- table(colData(sce)$sample_id,
             colData(sce)$cluster_id)
  
  metadata_ls <- list()
  
  for (i in 1:length(counts_list)) {
    df <- data.frame(cluster_sample_id = colnames(counts_list[[i]]))
    df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    idx <- which(colnames(t) == unique(df$cluster_id))
    cell_counts <- t[, idx]
    cell_counts <- cell_counts[cell_counts > 0]
    sample_order <- match(df$sample_id, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    df$cell_count <- cell_counts
    df <- plyr::join(df, metadata, 
                     by = intersect(names(df), names(metadata)))
    rownames(df) <- df$cluster_sample_id
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)
  }
  return(metadata_ls)
}

## deseq2_dea
# A wrapper function which does DEA using DESeq2 for pseudo-bulked single cell data. It automatically saves both significant and all DEGs in a 'pseudobulk' directory at the specified path. This code is originally form the HBC training guide, but was heavily adapted.
deseq2_dea <- function(cell_types, counts_ls, metadata_ls, group_oi, B, padj_cutoff = 0.05, path, shrinkage) {
  cell_type <- cell_types[1]
  print(cell_type)
  ifelse(!dir.exists(here(paste0(path, "03_dea/"))),
         dir.create(here(paste0(path, "03_dea/"))),
         print("Info: 03_dea directory already exists"))
  idx <- which(names(counts_ls) == cell_type)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ timepoint + group_id + timepoint:group_id)
  dds$group <- factor(paste0(dds$timepoint, dds$group_id))
  dds$group <- relevel(dds$group, ref = B)
  design(dds) <- ~ group
  dds <- DESeq(dds, quiet = TRUE)
  print(resultsNames(dds))
  name <- paste(c("group", group_oi, "vs", B), collapse = "_")
  print(name)
  res <- results(dds, name = name, alpha = 0.05)
  res <- lfcShrink(dds, coef = name, res = res, type = shrinkage, quiet = TRUE)
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  # save all results
  write.csv(res_tbl,
            here(paste0(path, "03_dea/", cell_type, "_", name, "_all_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
  
  # save sig results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            here(paste0(path, "03_dea/", cell_type, "_", name, "_sig_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
}

## make_names_valid
# Making names of interested columns valid (no spacing for example)
# Code adapted from MultiNicheNet
make_names_valid <- function(object) {
  SummarizedExperiment::colData(object)$ident <-
    SummarizedExperiment::colData(object)$ident %>%
    make.names()
  SummarizedExperiment::colData(object)$orig.ident <-
    SummarizedExperiment::colData(object)$orig.ident %>%
    make.names()
  SummarizedExperiment::colData(object)$sample <-
    SummarizedExperiment::colData(object)$sample %>%
    make.names()
  SummarizedExperiment::colData(object)$group <-
    SummarizedExperiment::colData(object)$group %>%
    make.names()
  return(object)
}

## multinichenet_wrapper
# A wrapper for original MultiNicheNet code
# Requires a single cell experiment object as input and returns a list of multinichenet outputs, which can be used for visualization. 
multinichenet_wrapper <- function(object, results_path, celltype_id, sample_id, group_id, lr_network, batches, contrasts_oi, contrast_tbl, covariates, empirical_pval, p_val_adj, cores_system, ligand_target_matrix, prioritizing_weights) {
  # get sender and receiver cell types from object ---------------
  senders_oi <- colData(object)[,celltype_id] %>% unique()
  receivers_oi <- colData(object)[,celltype_id] %>% unique()
  object <- object[, colData(object)[,celltype_id] %in% c(senders_oi, receivers_oi)]
  print("Grabbed senders and receivers")
  # determine whether all cell types have enough cells -----------
  min_cells <- 10
  abundance_expression_info <- get_abundance_expression_info(sce = object,
                                                             sample_id = sample_id,
                                                             group_id = group_id,
                                                             celltype_id = celltype_id,
                                                             min_cells = min_cells,
                                                             senders_oi = senders_oi,
                                                             receivers_oi = receivers_oi,
                                                             lr_network = lr_network,
                                                             batches = batches)
  plot1 <- abundance_expression_info$abund_plot_sample
  plot2 <- abundance_expression_info$abund_plot_group
  pdf(here(paste0(results_path, "/abundance_expression.pdf")))
  print(plot1)
  print(plot2)
  dev.off()
  print("Plotted abundance expression")
  # Get differential expression information
  DE_info <- get_DE_info(sce = object,
                         sample_id = sample_id,
                         group_id = group_id,
                         celltype_id = celltype_id,
                         batches = batches,
                         covariates = covariates,
                         contrasts_oi = contrasts_oi,
                         min_cells = min_cells)
  print("Calculated DE")
  # Calculate p-values
  if(empirical_pval == FALSE) {
    celltype_de <- DE_info$celltype_de$de_output_tidy
  } else {
    celltype_de <- DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>%
      dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
  }
  sender_receiver_de = combine_sender_receiver_de(
    sender_de = celltype_de,
    receiver_de = celltype_de,
    senders_oi = senders_oi,
    receivers_oi = receivers_oi,
    lr_network = lr_network)
  print("Calculated p-values")
  # Identify ligand activities
  fraction_cutoff <- 0.05
  n.cores <- min(cores_system, sender_receiver_de$receiver %>% unique() %>% length())
  if(length(unique(sender_receiver_de$receiver)) > n.cores) {
    print(paste0("Core req not met. Minimum number of cores needed: ",
                 length(unique(sender_receiver_de$receiver))))
    stop()
  } else {
    print("Calculating ligand activities")
    ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
      get_ligand_activities_targets_DEgenes(
        receiver_de = celltype_de,
        receivers_oi = receivers_oi,
        ligand_target_matrix = ligand_target_matrix,
        logFC_threshold = 0.50,
        p_val_threshold = 0.05,
        p_val_adj = p_val_adj,
        top_n_target = 250,
        verbose = FALSE, 
        n.cores = n.cores
      )))
  }
  # Prepare for prioritization
  print("Prioritizing interactions")
  sender_receiver_tbl <- sender_receiver_de %>%
    dplyr::distinct(sender, receiver)
  metadata_combined <- colData(object) %>% tibble::as_tibble()
  if(!is.na(batches)){
    grouping_tbl <- metadata_combined[,c(sample_id, group_id, batches)] %>%
      tibble::as_tibble() %>%
      dplyr::distinct()
    colnames(grouping_tbl) <- c("sample", "group", batches)
  } else {
    grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>%
      tibble::as_tibble() %>%
      dplyr::distinct()
    colnames(grouping_tbl) <- c("sample", "group")
  }
  # Prioritize interactions
  prioritization_tables <- suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    prioritizing_weights = prioritizing_weights,
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender
  ))
  # Prioritize by correlation coefficients
  print("Prior knowledge correlations")
  lr_target_prior_cor <- lr_target_prior_cor_inference(
    prioritization_tables$group_prioritization_tbl$receiver %>%
      unique(),
    abundance_expression_info,
    celltype_de,
    grouping_tbl,
    prioritization_tables,
    ligand_target_matrix,
    logFC_threshold = 0.50,
    p_val_threshold = 0.05,
    p_val_adj = p_val_adj)
  # Create combined output
  print("Creating multinichenet output")
  multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  )
  return(multinichenet_output)
}

## filter_nichenet
# A function that will filter a multinichenet output by pearson and spearman correlation values above/below 0.33/-0.33
filter_nichenet <- function(object) {
  # Filter correlated object by pearson and spearman correlations
  lr_target_prior_cor_filtered <- object$lr_target_prior_cor %>%
    inner_join(object$ligand_activities_targets_DEgenes$ligand_activities %>%
                 distinct(ligand, target, direction_regulation, contrast)) %>%
    inner_join(contrast_tbl) %>%
    filter(group == group_oi, receiver %in% receiver_oi, sender %in% sender_oi)
  
  lr_target_prior_cor_filtered_up <- lr_target_prior_cor_filtered %>%
    filter(direction_regulation == "up") %>%
    filter((rank_of_target < top_n_target) & (pearson > 0.33 | spearman > 0.33))
  
  lr_target_prior_cor_filtered_down <- lr_target_prior_cor_filtered %>%
    filter(direction_regulation == "down") %>%
    filter((rank_of_target < top_n_target) & (pearson < -0.33 | spearman < -0.33))
  
  lr_target_prior_cor_filtered <- bind_rows(lr_target_prior_cor_filtered_up,
                                            lr_target_prior_cor_filtered_down)
  # Create new column with ligand-receptor-target only info
  lr_target_prior_cor_filtered$ligand_receptor_target <- gsub(
    "(_[^_]+_)(.*?)(_[^_]+_)",
    "\\1\\4",
    lr_target_prior_cor_filtered$id_target,
    perl = TRUE)
  # Create new column with ligand-receptor only info
  lr_target_prior_cor_filtered$ligand_receptor <- gsub(
    "^([^_]*_[^_]*).*",
    "\\1",
    lr_target_prior_cor_filtered$id,
    perl = TRUE)
  return(lr_target_prior_cor_filtered)
}
