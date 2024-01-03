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






### LW code breakup
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

conditions <- c("6m_AD", "12m_AD", "6m_WT", "12m_WT")

condition_list <- object_list[grepl(conditions, names(object_list))]
make_seurat_obj <-   condition_list <- object_list[grepl("12m_AD", names(object_list))]
for (i in names(condition_list)) {
  sample_name <- basename(i)
  sample_name <- gsub("[[:punct:]]", "", sample_name)
  condition_list[[i]] <- RenameCells(condition_list[[i]],
                              add.cell.id = sample_name)
}
diseased_12m <- Merge_Seurat_List(condition_list)

