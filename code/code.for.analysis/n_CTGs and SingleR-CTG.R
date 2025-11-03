library(Seurat)


CTGs <- readRDS("../data/CTGs.rds")
seu  ## a seurat object



## n_CTGs ----------------------------------------------------------------------
genes  <- intersect(rownames(seu), CTGs) 
n.CTGs <- colSums(seu[["RNA"]]@data[genes, , drop = FALSE] > 0)  ## scRNA-seq
n.CTGs <- colSums(seu[["Spatial"]]@counts[genes, , drop = FALSE] > 0)  ## ST



## SingleR-CTG -----------------------------------------------------------------
top             <- n.CTGs >= quantile(n.CTGs, 0.95)
tail            <- n.CTGs <= quantile(n.CTGs, 0.05)
top.expr        <- seu[["RNA"]]@counts[, top, drop = FALSE] |> rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
tail.expr       <- seu[["RNA"]]@counts[, tail, drop = FALSE] |> rowSums() |> (\(x) log2(1e6 * x / sum(x) + 1))()
ref             <- cbind(top.expr, tail.expr)
pred            <- SingleR::SingleR(test = seu[["RNA"]]@data, ref = ref, labels = c("high", "low"))
malignant.score <- pred$scores[, "high"]