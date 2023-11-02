### Functions for CCC in 3xTgAD project 
# Tabea M. Soelter 

## remove_ambientRNA
# A function which removes ambient RNA from h5 files outputted from Cell Ranger for single cell data.
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
    filt_matrix <- Read10X_h5(paste0(i, "/filtered_feature_bc_matrix.h5"))
    raw_matrix <- Read10X_h5(paste0(i, "/raw_feature_bc_matrix.h5"))
    # create seurat object
    print("Making seurat object")
    object <- CreateSeuratObject(counts = filt_matrix)
    # make soup channel object
    print("Making soup channel object")
    sco <- SoupChannel(raw_matrix, filt_matrix)
    # get cluster info
    print("Get cluster info")
    object <- SCTransform(object, verbose = FALSE)
    object <- RunPCA(object, approx = FALSE, verbose = FALSE)
    object <- RunUMAP(object, dims = 1:30, verbose = FALSE)
    object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
    object <- FindClusters(object, verbose = FALSE)
    # ading metadata to soup channel object
    meta <- object@meta.data
    umap <- object@reductions$umap@cell.embeddings
    sco <- setClusters(sco, setNames(meta$seurat_clusters, rownames(meta)))
    sco <- setDR(sco, umap)
    # Analyzing the soup
    print("Profiling the soup")
    sco <- autoEstCont(sco)
    # Create integer matrix
    adjusted_matrix <- adjustCounts(sco, roundToInt = TRUE)
    # save
    print("Saving filtered objects")
    sample_name <- basename(i)
    DropletUtils::write10xCounts(paste0(outputs, sample_name), adjusted_matrix) 
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

## convert_human_gene_list
# Source: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/ 
# I adapted the function to use an archived version of ensembl, as there were mirror issues due to update in Feb 2023
# A function to convert human to mouse gene names 
convert_human_gene_list <- function(x) {
  require("biomaRt")
  human = useMart("ensembl",
                  dataset = "hsapiens_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl",
                  dataset = "mmusculus_gene_ensembl",
                  host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"),
                   filters = "hgnc_symbol",
                   values = x ,
                   mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}


## cell_cycle_effects
# A function which calculates and plots the effect of cell cycle on the data using a filtered seurat object as input. It also performs log normalization, scaling, and dimension reduction using PCA
cell_cylce_effects <- function(filtered_seurat){
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # convert human cell cylce markers to mouse -----
  s.genes <- convert_human_gene_list(cc.genes.updated.2019$s.genes)
  g2m.genes <- convert_human_gene_list(cc.genes.updated.2019$g2m.genes)
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m.genes,
                                      s.features = s.genes)
  filtered_seurat <- FindVariableFeatures(filtered_seurat,
                                          selection.method = "vst",
                                          verbose = FALSE)
  # scale data -----
  filtered_seurat <- ScaleData(filtered_seurat)
  # run pca -----
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
  seurat_object <- RunHarmony(seurat_object,
                              group.by.vars = "sample",
                              reduction = "pca",
                              dims.use = dims, assay.use = "RNA")
  # Here we use all PCs computed from Harmony for UMAP calculation -----
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
    dplyr::select(group_id, sample_id)
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






