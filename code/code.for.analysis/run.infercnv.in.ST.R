.libPaths("/home/fuchunyang/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(tidyverse)
library(Seurat)
library(infercnv)
library(foreach)
library(doParallel)



########## 1. run infercnv ########## 
load("./Proj_CTAs/input/ST/ST.data.RData")
df <- ST.data |> filter(calcu_CNV == "yes")


foreach(i = 1:nrow(df)) %do% {
  path   <- df$processed_path[i]
  sample <- df$unique_sample_id[i]
  
  seu           <- readRDS(path)
  out.dir       <- file.path("./Proj_CTAs/output/6. n_CTGs/ST/infercnv", sample)
  NormalCluster <- split(seu$immune_score, seu$seurat_clusters) |> sapply(mean) |> sort(decreasing = TRUE) |> head(n = 1) |> names()
  dir.create(out.dir)
  
  anno.file    <- seu[[]] |> select(seurat_clusters) |> rownames_to_column("barcode") 
  write.table(anno.file, file = paste(out.dir, "annotations_file.txt", sep = "/"), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix       = seu[["Spatial"]]@counts,
                                       gene_order_file         = "/home/fuchunyang/practice/inferCNV/input/hg38_gencode_v27.txt",
                                       annotations_file        = paste(out.dir, "annotations_file.txt", sep = "/"),
                                       ref_group_names         = NormalCluster,
                                       min_max_counts_per_cell = c(0, +Inf))  ## still remove cells with zero counts
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff                            = 0.1,
                                out_dir                           = out.dir,
                                cluster_by_groups                 = FALSE,
                                analysis_mode                     = "subclusters",
                                denoise                           = TRUE,
                                HMM                               = TRUE,
                                tumor_subcluster_partition_method = "random_trees",
                                HMM_type                          = "i6",
                                BayesMaxPNormal                   = 0,
                                num_threads                       = 30)
}



########## 2. CNV score ########## 
cl <- makeCluster(40)
registerDoParallel(cl)


load("./Proj_CTAs/input/ST/ST.data.RData")
path         <- "./Proj_CTAs/output/6. n_CTGs/ST/infercnv"
infercnv.df  <- foreach(i = dir(path), .combine = 'rbind') %do% {
  files             <- dir(file.path(path, i))
  hmm_finished      <- "07_tumor_subclustersHMMi6.rand_trees.random_trees.infercnv_obj" %in% files
  infercnv_finished <- "infercnv_obj.rds" %in% files
  data.frame(unique_sample_id = i, hmm_finished, infercnv_finished, infercnv_path = file.path(path, i))
}
infercnv.df  <- infercnv.df |> filter(infercnv_finished) |> left_join(ST.data, by = "unique_sample_id")


CNV.score.ls        <- foreach(i = 1:nrow(infercnv.df), .packages = c('infercnv', 'tibble', 'dplyr', 'foreach')) %dopar% {
  seu.path      <- infercnv.df$processed_path[i]
  infercnv.path <- infercnv.df$infercnv_path[i]
  seu           <- readRDS(seu.path)
  
  HMM.infercnv.obj       <- readRDS(file.path(infercnv.path, "17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.infercnv_obj"))
  cnv.score              <- data.frame(cnv_score = colSums(abs(HMM.infercnv.obj@expr.data - 3))) |> 
    rownames_to_column("barcode") |>
    mutate(normalized_cnv_score = cnv_score / nrow(HMM.infercnv.obj@expr.data))
  
  obs.grouping.ls        <- HMM.infercnv.obj@tumor_subclusters[["subclusters"]][["all_observations"]]
  obs.grouping.ls        <- lapply(obs.grouping.ls, names)
  names(obs.grouping.ls) <- as.character(seq_along(obs.grouping.ls))
  cnv.cluster            <- foreach(i = names(obs.grouping.ls), .combine = 'rbind') %do% {data.frame(barcode = obs.grouping.ls[[i]], cnv_cluster = i)}
  
  cnv.score <- cnv.score |> left_join(cnv.cluster) |> mutate(cnv_cluster = ifelse(is.na(cnv_cluster), "normal", cnv_cluster))
  stopifnot(all(cnv.score$barcode %in% colnames(seu)))
  cnv.score <- seu[[]] |> rownames_to_column("barcode") |> left_join(cnv.score)
}
names(CNV.score.ls) <- infercnv.df$unique_sample_id



