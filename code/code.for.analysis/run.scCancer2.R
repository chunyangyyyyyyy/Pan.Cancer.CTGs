.libPaths("/home/fuchunyang/R/x86_64-pc-linux-gnu-library/4.2/")  ## necessary for finding the packages used in the script.
library(Seurat)
library(foreach)
library(BiocParallel)


cl <- makeCluster(15)
registerDoParallel(cl)


load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
datasets <- dataset.info$dataset |> setdiff("Data_Brain_Richards2021_normalized")
param    <- MulticoreParam(workers = 20)


scCancer2.score.ls        <- bplapply(datasets, \(i) {
  tryCatch(
    expr = {
      prefix <- "/home/fuchunyang/Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object (normal included)/"
      file   <- paste0(i, ".rds")
      seu    <- readRDS(paste0(prefix, file))
      seu    <- seu |> FindVariableFeatures(nfeatures = 5000) |> ScaleData()
      
      load("/home/fuchunyang/practice/scCancer2/model.RData")
      source("/home/fuchunyang/practice/scCancer2/align_XGBoost.R")
      testdata <- t(seu[["RNA"]]@scale.data)
      testdata <- align_XGBoost(test = testdata, features = features)
      testdata <- xgboost::xgb.DMatrix(testdata)
      pred     <- predict(model.ref, testdata)
      df       <- data.frame(barcode = colnames(seu), sample = seu$sample, scCancer2_malignant_score = pred, row.names = NULL)
    },
    
    error = function(e) {e$message}
  )
}, BPPARAM = param)
names(scCancer2.score.ls) <- datasets
scCancer2.score.ls        <- foreach(i = datasets, .combine = 'c') %do% {
  scCancer2.score <- scCancer2.score.ls[[i]]
  ls              <- split(scCancer2.score[, c("barcode", "scCancer2_malignant_score")], 
                           scCancer2.score$sample)
  names(ls)       <- paste(i, names(ls), sep = "_")
  ls
}
bpstop(param)