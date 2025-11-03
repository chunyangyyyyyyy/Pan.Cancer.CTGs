.libPaths("/home/fuchunyang/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(Seurat)
library(tidyverse)
library(foreach)
library(doParallel)



cl <- makeCluster(20)
registerDoParallel(cl)


load("/home/fuchunyang/Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
datasets <- dataset.info$dataset |> setdiff("Data_Brain_Richards2021_normalized")


## 1. run ----------------------------------------------------------------------
foreach(i = datasets, .packages = 'Seurat') %dopar% {
  prefix <- "/home/fuchunyang/Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object (normal included)/"
  file   <- paste0(i, ".rds")
  seu    <- readRDS(paste0(prefix, file))
  
  seu.ls        <- SplitObject(seu, split.by = "sample")
  names(seu.ls) <- paste(i, names(seu.ls), sep = "_")
  
  lapply(names(seu.ls), \(ds) {
    seu    <- seu.ls[[ds]]
    counts <- seu[["RNA"]]@counts
    path   <- "/home/fuchunyang/Proj_CTAs/output/6. n_CTGs/scRNA-seq/ikarus"
    path   <- file.path(path, ds)
    dir.create(path)
    setwd(path)
    write.table(t(as.matrix(counts)), file = "expr.txt", sep = "\t", quote = FALSE, col.names = NA)
    write.table("Fished!", file = "finish.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
    gc()
    
    tryCatch({
      system("/home/fuchunyang/practice/ikarus/run_ikarus.sh")
    }, error = function(e) {
      write(paste("Error in directory: ", path, " - ", e$message), file = "error.log.txt", append = TRUE)
    })
  })
}



## 2. score --------------------------------------------------------------------
prefix <- "/home/fuchunyang/Proj_CTAs/output/6. n_CTGs/scRNA-seq/ikarus"
df     <- data.frame(dataset_sample = dir(prefix), path = dir(prefix, full.names = TRUE))
df     <- df |> mutate(finished = sapply(path, \(i) "prediction.csv" %in% dir(file.path(i, "out", "predict"))))


all.ds <- sample.info |> mutate(dataset_sample = paste(dataset, sample, sep =  "_")) |> filter(dataset != "Data_Brain_Richards2021_normalized")
setdiff(all.ds$dataset_sample, df$dataset_sample)  ## all 95 datasets included
df     <- df |> 
  filter(finished) |> 
  mutate(pred_path = file.path(path, "out", "predict", "prediction.csv"),
         expr_path = file.path(path, "expr.txt")) |> 
  left_join(all.ds, by = "dataset_sample")


ikarus.score.ls        <- foreach(i = 1:nrow(df)) %dopar% {
  expr <- data.table::fread(df$expr_path[i])
  pred <- read.csv(df$pred_path[i])
  stopifnot(nrow(expr) == nrow(pred))
  data.frame(barcode = expr$V1, ikarus_malignant_score = pred$final_pred_proba_Tumor)
}
names(ikarus.score.ls) <- df$dataset_sample