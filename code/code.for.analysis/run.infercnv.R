.libPaths("/home/fuchunyang/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(foreach)
library(doParallel)
library(tidyverse)
library(Seurat)
library(BiocParallel)
library(infercnv)



load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
coarse.cell.type.ls <- readRDS("./Proj_CTAs/output/2. cell annotation/coarse.cell.type.ls.rds")


########## 1. run: coarse cell type (91) ########## 
df <- foreach(i = names(coarse.cell.type.ls), .combine = 'rbind') %do% {
  coarse.cell.type.ls[[i]] |> 
    filter(my_coarse_cell_type != "Unknown") |>
    group_by(sample) |> 
    count(my_coarse_cell_type) |>
    mutate(dataset = i, sample = as.character(sample), dataset_sample = paste(dataset, sample, sep = "_"))
}
df <- df |> left_join(sample.info, by = c("dataset", "sample"))


## both tumor and normal samples
samples.contain.immume.or.stromal       <- df |> filter(my_coarse_cell_type %in% c("Stromal", "Immune")) |> pull(dataset_sample) |> unique()
samples.contain.more.than.one.cell.type <- df |> group_by(dataset_sample) |> count() |> filter(n > 1) |> pull(dataset_sample) |> unique()
samples.kept                            <- intersect(samples.contain.immume.or.stromal, samples.contain.more.than.one.cell.type)


args.df1 <- df |> filter(dataset_sample %in% samples.kept, my_coarse_cell_type %in% c("Stromal", "Immune")) |> group_by(dataset_sample) |> slice_max(order_by = n, n = 1, with_ties = FALSE)
args.df1 <- args.df1 |> ungroup() |> filter(n >= 10) |> left_join(dataset.info, by = "dataset") |> mutate(cutoff = ifelse(tech == "SmartSeq2", 1, 0.1))


param   <- MulticoreParam(workers = 20)
results <- bplapply(unique(args.df1$dataset), function(dataset) {
  prefix           <- "./Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object (normal included)/"
  file             <- paste0(dataset, ".rds")
  seu              <- readRDS(paste0(prefix, file))
  coarse.cell.type <- coarse.cell.type.ls[[dataset]]
  stopifnot(all(rownames(coarse.cell.type) == colnames(seu)))
  seu$cell_type    <- coarse.cell.type$my_coarse_cell_type
  seu.ls           <- SplitObject(seu, split.by = "sample")
  
  
  samples         <- args.df1$sample[args.df1$dataset == dataset]
  refs            <- args.df1$my_coarse_cell_type[args.df1$dataset == dataset]
  cutoffs         <- args.df1$cutoff[args.df1$dataset == dataset]
  dataset.samples <- args.df1$dataset_sample[args.df1$dataset == dataset]
  
  
  foreach(s = samples, r = refs, c = cutoffs, ds = dataset.samples) %do% {
    seu     <- seu.ls[[s]]
    out.dir <- file.path("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/infercnv", ds)
    dir.create(out.dir)
    
    anno.file    <- seu[[]] |> select(cell_type) |> rownames_to_column("barcode") 
    write.table(anno.file, file = paste(out.dir, "annotations_file.txt", sep = "/"), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix       = seu[["RNA"]]@counts,
                                         gene_order_file         = "/home/fuchunyang/practice/inferCNV/input/hg38_gencode_v27.txt",
                                         annotations_file        = paste(out.dir, "annotations_file.txt", sep = "/"),
                                         ref_group_names         = r,
                                         min_max_counts_per_cell = c(0, +Inf))  ## avoid filter out cells of some datasets
    
    ## only for calculating cnv score
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff            = c,
                                  out_dir           = out.dir,
                                  cluster_by_groups = FALSE,
                                  denoise           = FALSE,
                                  HMM               = FALSE,
                                  up_to_step        = 14,
                                  num_threads       = 1,
                                  no_plot           = TRUE,
                                  write_expr_matrix = TRUE)
    
    saveRDS(infercnv_obj, paste(out.dir, "infercnv_obj.rds", sep = "/"))
  }
  
  gc()
}, BPPARAM = param)
bpstop(param)



########## 2. run: malignant (yes/no) (4) ########## 
## only tumor samples
datasets <- c("Data_Lung_Maynard2020_LUAD", "Data_Lung_Maynard2020_LUSC", "Data_Skin_AlvarezBreckenridge2022_log2TP100K", "Data_Neuroendocrine_Dong2020_Group2")
args.df2 <- tumor.samples |> 
  filter(dataset %in% datasets) |> 
  left_join(dataset.info, by = "dataset") |>
  mutate(dataset_sample = paste(dataset, sample, sep = "_"), cutoff = ifelse(tech == "SmartSeq2", 1, 0.1))


param   <- MulticoreParam(workers = 4)
results <- bplapply(datasets, function(dataset) {
  prefix <- "./Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object (normal included)/"
  file   <- paste0(dataset, ".rds")
  seu    <- readRDS(paste0(prefix, file))
  seu.ls <- SplitObject(seu, split.by = "sample")
  
  samples         <- args.df2$sample[args.df2$dataset == dataset]
  cutoffs         <- args.df2$cutoff[args.df2$dataset == dataset]
  dataset.samples <- args.df2$dataset_sample[args.df2$dataset == dataset]
  
  foreach(s = samples, c = cutoffs, ds = dataset.samples) %do% {
    seu     <- seu.ls[[s]]
    out.dir <- file.path("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/infercnv", ds)
    dir.create(out.dir)
    
    anno.file    <- seu[[]] |> select(malignant) |> rownames_to_column("barcode") 
    write.table(anno.file, file = paste(out.dir, "annotations_file.txt", sep = "/"), sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix       = seu[["RNA"]]@counts,
                                         gene_order_file         = "/home/fuchunyang/practice/inferCNV/input/hg38_gencode_v27.txt",
                                         annotations_file        = paste(out.dir, "annotations_file.txt", sep = "/"),
                                         ref_group_names         = "no",
                                         min_max_counts_per_cell = c(0, +Inf))  ## avoid filter out cells of some datasets
    
    ## only for calculating cnv score
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff            = c,
                                  out_dir           = out.dir,
                                  cluster_by_groups = FALSE,
                                  denoise           = FALSE,
                                  HMM               = FALSE,
                                  up_to_step        = 14,
                                  num_threads       = 1,
                                  no_plot           = TRUE,
                                  write_expr_matrix = TRUE)
    
    saveRDS(infercnv_obj, paste(out.dir, "infercnv_obj.rds", sep = "/"))
  }
  
  gc()
}, BPPARAM = param)
bpstop(param)


