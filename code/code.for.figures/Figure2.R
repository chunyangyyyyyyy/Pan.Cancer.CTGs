library(tidyverse)
library(patchwork)
library(Seurat)
library(ggrepel)
library(foreach)
library(ComplexHeatmap)


source("./Proj_CTAs/code/visualization.R")
load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")


calculate.tau <- function(x) {
  stopifnot(all(x >= 0))
  normalized.x <- x / max(x)
  tau          <- sum(1 - normalized.x) / (length(x) - 1)
  return(tau)
}


## Fig. 2A ---------------------------------------------------------------------
# ## TNBC and non-TNBC samples in TCGA
# BRCA.info        <- data.table::fread("./Datas/UCSC_Xena/TCGA_TARGET_GTEx/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix")
# BRCA.info        <- BRCA.info |>
#   select(sampleID, ER = breast_carcinoma_estrogen_receptor_status, PR = breast_carcinoma_progesterone_receptor_status, HER2 = lab_proc_her2_neu_immunohistochemistry_receptor_status) |>
#   filter(ER %in% c("Negative", "Positive"), PR %in% c("Negative", "Positive"), HER2 %in% c("Negative", "Positive"),
#          as.numeric(str_sub(sampleID, 14)) < 10)  ## only keep tumour samples
# TNBC.samples     <- BRCA.info |> filter(ER == "Negative", PR == "Negative", HER2 == "Negative") |> pull(sampleID)
# non.TNBC.samples <- setdiff(BRCA.info$sampleID, TNBC.samples)
# 
# 
# load("./Datas/UCSC_Xena/tcga.batch.normalized.RData")
# coldata <- tcga.sample.info |> mutate(group = case_match(sample, TNBC.samples ~ "TNBC", non.TNBC.samples ~ "non-TNBC", .default = "Other Cancers")) |> column_to_rownames("sample")
# coldata <- coldata |> mutate(group = factor(group, levels = c("TNBC", "non-TNBC", "Other Cancers")))
# seu     <- CreateSeuratObject(counts = tcga.log2.batch.normalized.exp, meta.data = coldata)
# 
# 
# CTGs <- intersect(rownames(seu), pan.sample.CTGs)
# seu  <- ScaleData(seu, features = CTGs)
# seu  <- RunPCA(seu, features = CTGs)
# seu  <- RunUMAP(seu, dims = 1:30, reduction = "pca")



load("../data/fig.2a.RData")  ## seu
cancer.info <- TCGAbiolinks::getGDCprojects() |> filter(grepl("^TCGA", id)) |> select(id, cancer = name, abbr = tumor)
cancers     <- cancer.info$abbr
colors      <- setNames(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Set3"))(length(cancers)), cancers)
p1          <- DimPlot(seu, group.by = "abbr", label = TRUE, repel = TRUE, pt.size = 0.01, raster = FALSE, label.size = 6.5) + 
  scale_color_manual(values = colors) + 
  labs(title = NULL) +
  tidydr::theme_dr() +
  coord_fixed() +
  theme(aspect.ratio = 1, legend.position = "bottom")
p2          <- DimPlot(seu, group.by = "group", pt.size = 0.01, raster = FALSE) + 
  scale_color_manual(values = c("TNBC" = "skyblue", "non-TNBC" = "orange", "Other Cancers" = "lightgrey")) + 
  labs(title = NULL) +
  theme_void() +
  NoLegend() +
  coord_fixed() +
  theme(aspect.ratio = 1, panel.border = element_rect(linetype = "dotted", linewidth = 0.5, fill = NA)) 
label.df    <- p2[[1]][["data"]] |> group_by(group) |> summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)) |> filter(group != "Other Cancers")
p2          <- p2 + 
  geom_text_repel(data = label.df, aes(UMAP_1, UMAP_2, label = group), 
                  size = 4.5, point.padding = 0.5, box.padding = 0.6, min.segment.length = 0)


p1 + inset_element(p2, left = 0.65, bottom = 0.6, right = 1, top = 1) & 
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 4)))



## Fig. 2B ---------------------------------------------------------------------
# ccle.log2.tpm        <- data.table::fread("./Datas/CCLE/DepMap_22Q2/CCLE_expression.csv") |> column_to_rownames("V1")
# names(ccle.log2.tpm) <- str_extract(names(ccle.log2.tpm), ".+(?=\\s)")
# ccle.info            <- read.csv("./Datas/CCLE/DepMap_22Q2/sample_info.csv")


load("../fig.2b.RData")
bc         <- ccle.info |> filter(primary_disease == "Breast Cancer", lineage_molecular_subtype != "")
bc$subtype <- factor(bc$lineage_molecular_subtype, levels = c("basal", "basal_A", "basal_B", "luminal", "HER2_amp", "luminal_HER2_amp"))


plots <- foreach(i = c("ELL3", "TEKT5", "ROPN1")) %do% {
  cbind(bc, tpm = ccle.log2.tpm[bc$DepMap_ID, i]) |>
    ggplot(aes(x = subtype, y = tpm)) +
    geom_boxplot(outliers = FALSE, color = "#2F4F4F", fill = "#F0F8FF", linewidth = 0.3, alpha = 0.7) +
    geom_jitter(size = 0.01, width = 0.1, color = "#4682B4") +
    labs(x = NULL, y = expression(Log[2](TPM + 1)), title = i) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"), axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
}
plots |> wrap_plots(ncol = 3)



## Fig. 2C-D -------------------------------------------------------------------
CTG.class        <- unlist(pan.sample.CTGs.ls) |> table() |> as.data.frame()
names(CTG.class) <- c("CTG", "n")
CTG.class        <- CTG.class |> 
  mutate(CTG   = as.character(CTG),
         class = case_when(
           n >= 11 ~ "pan-cancer",
           n <= 3  ~ "cancer-specific",
           .default = "other"
         ))


##### (1) C #####
df <- CTG.class |> 
  count(n, name = "count") |>
  mutate(x     = factor(as.character(n), levels = 1:13),
         class = case_when(
           n >= 11 ~ "pan-cancer",
           n <= 3  ~ "cancer-specific",
           .default = "other"
         ))
p  <- ggplot(df, aes(x = x, y = count, fill = class)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "black", linewidth = 0.2) +
  geom_text(aes(label = count), vjust = -0.5) +
  scale_y_continuous(expand =  expansion(mult = c(0, 0.05))) +
  labs(x = "Number of cancer types", y = "Frequency", fill = "CTG class") +
  theme_classic() +
  scale_fill_manual(values = c("cancer-specific" = "#5B9BD5", "pan-cancer" = "#ED7D31", "other" = "#E7E6E6")) +
  theme(legend.position = c(0.9, 0.8),
        legend.box.background = element_rect(color = "black", linewidth = 1, linetype = "solid"),
        legend.key = element_rect(color = "gray70"),
        legend.key.size = unit(0.4, "cm"))
p


##### (2) D #####
# load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# load("./Proj_CTAs/output/3. define CTGs/pct.summary.RData")
# 
# 
# mat             <- pct.mat[pan.sample.CTGs, ]
# mat[is.na(mat)] <- 0  ## replace NAs to zero
# 
# 
# col.df <- tumor.samples |> 
#   inner_join(dataset.info.10x) |> 
#   left_join(cancer.abbr) |> 
#   mutate(sample = paste(dataset, sample, sep = "_"), abbr = factor(abbr, levels = cancer.abbr$abbr)) |> 
#   column_to_rownames("sample")
# 
# 
# 
# median.mat                    <- foreach(i = cancer.abbr$cancer, .combine = 'rbind') %do% {data.frame(abbr = cancer.abbr$abbr[cancer.abbr$cancer == i], CTG = pct.summary.each.cancer.10x[[i]]$CTG, value = pct.summary.each.cancer.10x[[i]]$pct_median)} |> pivot_wider(names_from = abbr, values_from = value) |> column_to_rownames("CTG")
# median.mat[is.na(median.mat)] <- 0
# row.df                        <- foreach(i = CTG.class$CTG[CTG.class$class == "cancer-specific"], .combine = 'rbind') %do% {
#   data.frame(CTG = i, abbr = names(which.max(median.mat[i, ])), pct = max(median.mat[i, ]))
# } |> 
#   mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) |>
#   arrange(abbr, desc(pct)) |>
#   column_to_rownames("CTG")
# 
# 
# mat <- mat[rownames(row.df), rownames(col.df)]


load("../data/fig.2d.RData")
selected.genes <- c("IL13RA2", "ROPN1", "PRAME", "C1orf94", "MDGA2", "PRDM7", "NLRP7", "TEKT5", "HTR7", "PLAC1", "CSAG1", "IGF2BP1", "MAGEA1", "TRIM71", "LRRIQ4")
index          <- match(selected.genes, rownames(row.df))
label.anno     <- rowAnnotation(mark = anno_mark(at = index, labels = selected.genes,link_gp = gpar(lwd = 1), labels_gp = gpar(fontsize = 7, fontface = "italic")))


scaled.mat <- t(scale(t(mat)))
quantile   <- quantile(scaled.mat, probs = c(0.01, 0.99), na.rm = TRUE)
col.fun    <- circlize::colorRamp2(seq(quantile[1], quantile[2], length.out = 100), colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlGn")))(100))
col.ha     <- HeatmapAnnotation(Cancer = col.df$abbr, col = list(Cancer = cancer.colours), show_legend = TRUE, show_annotation_name = TRUE)
row.ha     <- rowAnnotation(cancer = row.df$abbr, col = list(cancer = cancer.colours), show_legend = FALSE, show_annotation_name = FALSE)


p <- Heatmap(scaled.mat, col = col.fun, name = "z-score", 
             cluster_rows = FALSE, cluster_columns = FALSE,
             show_row_names = FALSE, show_column_names = FALSE,
             column_split = col.df$abbr,  row_split = row.df$abbr,
             top_annotation = col.ha, left_annotation = row.ha, right_annotation = label.anno,
             column_title = " ", column_title_gp = gpar(col = "black", fontface = "bold"),
             row_names_gp = gpar(col = row.df$color, fontsize = 8, fontface = "italic", lwd = 0.01),
             row_title = " ",
             border = TRUE, border_gp = gpar(col = "black", lwd = 0.01))
p <- draw(p)



## Fig. 2E ---------------------------------------------------------------------
# load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# CTGs.df <- CTGs.df |> 
#   mutate(is_CTG        = ifelse(symbol %in% pan.sample.CTGs, "yes", "no"),
#          chr_group     = case_when(
#            grepl("^X", location) ~ "X chromosome", 
#            grepl("^Y", location) ~ "Y chromosome", 
#            .default = "autosome"
#          )) 
# df      <- CTGs.df |> filter(chr_group != "Y chromosome")
# fig.2e  <- df


df    <- fig.2e
p.val <- fisher.test(table(df$chr_group, df$is_CTG))$p.value
p     <- df |>
  ggplot(aes(x = chr_group, fill = is_CTG)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_manual(values = c("yes" = "#FBB4AE", "no" = "#B3CDE3"), name = "is_CTG") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent_format(), expand =  expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Proportion of CTGs", title = NULL) + 
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  annotate("text", x = 1.5, y = 1.05, 
           size = 4, color = "black", lineheight = 1,
           label = paste0("Fisher's exact test, p = ", formatC(p.val, format = "g", digits = 2)))
p


## Fig. 2F ---------------------------------------------------------------------
# load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# load("./Proj_CTAs/output/3. define CTGs/pct.summary.RData")
# CTGs.df    <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# x.genes    <- CTGs.df$symbol[grepl("^X", CTGs.df$location)]
# auto.genes <- CTGs.df$symbol[grepl("^\\d", CTGs.df$location)]
# 
# 
# tau.df <- foreach(i = cancers, .combine = 'rbind') %do% {
#   tau <- pct.summary.each.cancer.10x[[i]] |> select(-pct_median) |> column_to_rownames("CTG") |> apply(1, calculate.tau)
#   data.frame(gene = names(tau), tau, cancer = i, row.names = NULL)
# }
# tau.df <- tau.df |> mutate(group = case_match(gene, x.genes ~ "X-CTGs (Identified)", auto.genes ~ "Autosomal CTGs (Identified)")) |> filter(!is.na(group)) |> left_join(cancer.abbr)
# fig.2f <- tau.df


tau.df <- fig.2f
p      <- tau.df |> filter(gene %in% pan.sample.CTGs) |>
  ggplot(aes(x = abbr, y = tau, color = group)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 1.2, alpha = 0.3) +
  scale_color_manual(values = c("X-CTGs (Identified)" = "#D53E4F", "Autosomal CTGs (Identified)" = "#4575B4"), name = "CTG Group") +
  scale_x_discrete(limit = cancer.abbr$abbr) +
  theme_bw() +
  labs(x = NULL, y = "Tau", title = "discovery cohort") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5), legend.position = "top") +
  ggpubr::stat_compare_means(label = "p.signif", hide.ns = TRUE, show.legend = FALSE)
p
