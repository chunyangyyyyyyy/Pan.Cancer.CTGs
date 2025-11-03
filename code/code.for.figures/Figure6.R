library(Seurat)
library(foreach)
library(tidyverse)


## Fig. 6A ---------------------------------------------------------------------
# source("./Proj_CTAs/code/visualization.R")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")  ## c("sample.ls", "tumor.sample.ls", "gene.info", "log2.tpm")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# CTGs.df <- CTGs.df |> filter(symbol %in% pan.sample.CTGs) |> rename(CTG = symbol) |> inner_join(gene.info, by = "ensembl")
# 
# 
# CTG.auc.bulk.ls        <- foreach(i = names(sample.ls)) %do% {
#   samples <- sample.ls[[i]]
#   y       <- factor(samples$group, levels = c("normal", "tumour"))
#   auc.df  <- foreach(g = CTGs.df$id, .combine = 'rbind') %do% {
#     x     <- log2.tpm[g, samples$sample]
#     p.val <- wilcox.test(x ~ y)$p.value
#     auc   <- pROC::roc(y, x, direction = "<")$auc
#     data.frame(id = g, auc, p.val)
#   }
#   auc.df  <- auc.df |> mutate(cancer = i) |> left_join(CTGs.df, by = "id")
#   auc.df  <- auc.df |> filter(CTG %in% pan.sample.CTGs.ls[[i]]) |> mutate(p.adj = p.adjust(p.val, method = "fdr"))
# }
# names(CTG.auc.bulk.ls) <- names(sample.ls)
# 
# 
# gtex.log2.tpm    <- data.table::fread("./Datas/UCSC_Xena/GTEx/gtex_RSEM_gene_tpm.gz") |> column_to_rownames("sample") |> as.matrix()
# gtex.tpm         <- 2^gtex.log2.tpm - 0.001
# all(CTGs.df$id %in% rownames(gtex.tpm))  ## TRUE
# gtex.sample.info <- data.table::fread("./Datas/UCSC_Xena/GTEx/GTEX_phenotype.gz")
# gtex.sample.info <- gtex.sample.info |> mutate(tissue = sapply(strsplit(`body_site_detail (SMTSD)`, " - "), `[[`, 1)) |> filter(Sample %in% colnames(gtex.tpm))
# gtex.sample.ls   <- split(gtex.sample.info$Sample, gtex.sample.info$tissue)
# gtex.median.tpm  <- sapply(gtex.sample.ls, \(i) {gtex.tpm[CTGs.df$id, i] |> apply(1, median)})
# gtex.median.tpm  <- gtex.median.tpm[, colnames(gtex.median.tpm) %in% c("Cells", "Testis") == FALSE]
# gtex.max.tpm     <- data.frame(CTG        = CTGs.df$CTG, 
#                                max.tpm    = apply(gtex.median.tpm, 1, \(x) sort(x, decreasing = TRUE)[1]), 
#                                second.tpm = apply(gtex.median.tpm, 1, \(x) sort(x, decreasing = TRUE)[2]))
# 
# 
# load("./Proj_CTAs/output/3. define CTGs/pct.summary.RData")
# cancers               <- names(pct.summary.each.cancer.10x)
# auc.tpm.pct.ls        <- foreach(i = cancers) %do% {
#   ls <- list(
#     CTG.auc.bulk.ls[[i]] |> select(cancer, CTG, ensembl, auc.bulk = auc, p.val.bulk = p.val, p.adj.bulk = p.adj),
#     gtex.max.tpm,
#     pct.summary.each.cancer.10x[[i]] |> select(CTG, pct = pct_median)
#   )
#   df <- Reduce(inner_join, ls)
# }
# names(auc.tpm.pct.ls) <- cancers
# 
# 
# dot.df <- foreach(i = auc.tpm.pct.ls, .combine = 'rbind') %do% {i |> filter(auc.bulk > 0.8, p.adj.bulk < 0.01, max.tpm < 1) |> select(CTG, cancer, auc.bulk, pct, max.tpm)}
# orders <- c("CKAP2L", "KIF18A", "BRIP1", "AUNIP", "HSF2BP", "IGF2BP3", "FBXO43", "NACA2",
#             "EIF5AL1", "RDM1", "STK31", "KDM4D", "SMC1B", "SPATA12", "CSAG1", "LRRIQ3", 
#             "C1orf94", 
#             "CT83", 
#             "TEKT5",
#             "LRRIQ4", "PIWIL1", 
#             "HTR7", "PLAC1",
#             "MYO3A",
#             "IGF2BP1", "MAGEA1", "TRIM71",
#             "ACRV1", "CCDC62", "EFCAB5", "HMX1", "PLAC8L1", "PPEF1",
#             "NLRP7", "PLSCR2",
#             "MDGA2", "PRDM7")
# dot.df <- dot.df |> 
#   left_join(cancer.abbr) |> 
#   left_join(CTGs.df) |>
#   mutate(abbr = factor(abbr, levels = cancer.abbr$abbr), CTG = factor(CTG, levels = orders))
# fig.6a <- dot.df


dot.df <- fig.6a
p      <- dot.df |>
  mutate(CTG = factor(CTG, levels = rev(levels(dot.df$CTG)))) |>
  ggplot(aes(x = abbr, y = CTG)) +
  geom_point(aes(size = pct, color = auc.bulk), shape = 19, alpha = 0.8) +
  scale_color_gradientn(name = "ROC-AUC", 
                        colors = c("#4477AA", "#77AADD", "#B4D8A9", "#FFDD77", "#FF9966"),
                        limits = c(0.8, 1),
                        guide = guide_colorbar(order = 2)) +
  scale_size_continuous(name = "Malignant expr. (%)",
                        guide = guide_legend(order = 1)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x.top = element_text(angle = 30, hjust = 0, vjust = 0),
        axis.text.y.right = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank()) +
  geom_hline(yintercept = 1:36 + 0.5, linewidth = 0.3, color = "grey90") +
  geom_vline(xintercept = 1:12 + 0.5, linewidth = 0.3, color = "grey90")


left.label <- dot.df |> distinct(CTG, max.tpm) |> 
  mutate(max.tpm = ifelse(max.tpm < 0, 0, max.tpm)) |>
  ggplot(aes(x = 1, y = CTG)) + 
  geom_bar(aes(fill = max.tpm, color = max.tpm), stat = "identity") +
  scale_fill_viridis_c(direction = -1, option = "A", limits = c(0, 1),
                       name = "Max non-testis normal expr. (TPM)") +
  scale_color_viridis_c(direction = -1, option = "A", limits = c(0, 1),
                        name = "Max non-testis normal expr. (TPM)") +
  theme_void() +
  theme(legend.position = "bottom")


p |> aplot::insert_left(left.label + theme(legend.position = "none"), width = 0.025)
