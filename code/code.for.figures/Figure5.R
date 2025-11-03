library(tidyverse)
library(foreach)
library(Seurat)
library(patchwork)
library(cowplot)


load("./Proj_CTAs/input/ST/ST.data.RData")
ST.data <- ST.data |> 
  mutate(cancer     = ifelse(cancer == "OSCC", "HNSCC", cancer),
         with_image = ifelse(folder %in% c("GSE230282_PRAD", "GSE200278_SKCM"), "no", "yes")) |>
  filter(cancer %in% c("GIST", "lymphoma", "OCS") == FALSE)



## Fig. 5A ---------------------------------------------------------------------
# df     <- foreach(i = grep("^GSE208253_OSCC", names(ST.meta.ls), value = TRUE), .combine = 'rbind') %do% {
#   ST.meta.ls[[i]] |> 
#     select(pathologist_anno.x, n_CTGs, ESTIMATE_purity) |> 
#     filter(!is.na(pathologist_anno.x)) |>
#     mutate(sample = str_sub(i, start = 11), y = ifelse(pathologist_anno.x == "SCC", 1, 0))
# } 
# df     <- df |> mutate(sample    = factor(sample, levels = paste0("OSCC_sample_", 1:12)),
#                        malignant = ifelse(pathologist_anno.x == "SCC", "yes", "no"))
# fig.5a <- df


df <- fig.5a
p  <- ggplot() + 
  gghalves::geom_half_violin(data = df |> filter(malignant == "no"), aes(x = sample, y = log10(n_CTGs + 1)), 
                             fill = "#377EB8", color = "#377EB8", side = "l", nudge = 0.01, alpha = 0.5) +
  gghalves::geom_half_violin(data = df |> filter(malignant == "yes"), aes(x = sample, y = log10(n_CTGs + 1)), 
                             fill = "#E41A1C", color = "#E41A1C", side = "r", nudge = 0.01, alpha = 0.5) +
  gghalves::geom_half_boxplot(data = filter(df, malignant == "no"), aes(x = sample, y = log10(n_CTGs + 1)), 
                              fill = "white", color = "grey40", side = "l", nudge = 0.03, width = 0.3, outlier.size = 0.1) +
  gghalves::geom_half_boxplot(data = filter(df, malignant == "yes"), aes(x = sample, y = log10(n_CTGs + 1)), 
                              fill = "white", color = "grey40", side = "r", nudge = 0.03, width = 0.3, outlier.size = 0.1) +
  ggpubr::stat_compare_means(data = df, aes(x = sample, y = log10(n_CTGs + 1), group = malignant), label = "p.signif", hide.ns = TRUE) +
  labs(x = NULL, y = expression(log[10]("n_CTGs"))) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 30))


legend.p <- ggplot(df, aes(x = malignant, y = n_CTGs, fill = malignant, color = malignant)) + 
  geom_violin(alpha = 0.5) +
  scale_color_manual(values = c("yes" = "#E41A1C", "no" = "#377EB8")) +
  scale_fill_manual(values = c("yes" = "#E41A1C", "no" = "#377EB8")) +
  theme_classic()
legend   <- cowplot::get_legend(legend.p)


p + wrap_elements(full = legend) + plot_layout(widths = c(9, 1))



## Fig. 5B ---------------------------------------------------------------------
# df            <- ST.data |> filter(folder == "Mendeley_skrx2fz79n_HCC")
# seu.ls        <- foreach(i = df$processed_path) %do% {readRDS(i)}
# names(seu.ls) <- df$unique_sample_id
# 
# 
# ## the annotation of the two samples were consistent with the original paper
# samples.kept <- c("Mendeley_skrx2fz79n_HCC_P1T", "Mendeley_skrx2fz79n_HCC_P3T")
# seu.ls.kept  <- seu.ls[samples.kept]
# cell.types   <- lapply(seu.ls.kept, \(x) as.character(unique(x$DefineTypes))) |> unlist() |> unique()
# cell.levels  <- c("Malignant hepatocyte", "Normal hepatocyte", "Immune/Fibroblast", "Myofibroblast/Pericyte", "SPP1+ Macrophage/CAF", "Cholangiocyte")
# seu.ls.kept  <- lapply(seu.ls.kept, \(x) {
#   x$cell_type <- case_match(
#     as.character(x$DefineTypes),
#     c("HMGB2 malignant hepatocyte", "MAP3K12 malignant hepatocyte", "Malignant hepatocyte", "malignant hepatocyte") ~ "Malignant hepatocyte",
#     "SPP1_Macrophage/CAF"    ~ "SPP1+ Macrophage/CAF",
#     "Myofibroblast_Pericyte" ~ "Myofibroblast/Pericyte",
#     "Immune_Fibroblast"      ~ "Immune/Fibroblast",
#     "hepatocyte"             ~ "Normal hepatocyte",
#     "Cholangiocyte"          ~ "Cholangiocyte")
#   x$cell_type <- factor(x$cell_type, levels = cell.levels)
#   x
# })


load("../data/fig.5b.RData")
plots <- foreach(i = names(seu.ls.kept)) %do% {
  normal.cells <- c("Normal hepatocyte", "Immune/Fibroblast", "Myofibroblast/Pericyte", "SPP1+ Macrophage/CAF", "Cholangiocyte")
  ls           <- lapply(intersect(normal.cells, seu.ls.kept[[i]]$cell_type), append, "Malignant hepatocyte")
  y.max        <- max(seu.ls.kept[[i]]$n_CTGs)
  yintercept   <- seu.ls.kept[[i]][[]] |> filter(cell_type == "Malignant hepatocyte") |> pull(n_CTGs) |> median()
  
  p <- VlnPlot(seu.ls.kept[[i]], features = "n_CTGs", group.by = "cell_type", y.max = y.max * (1 + length(ls) * 0.11), pt.size = 0) +
    geom_jitter(width = 0.08, color = "black", size = 0.1) +
    stat_summary(fun = median, geom = "point", size = 4, colour = "white", shape = 18) +
    labs(x = NULL, subtitle = i, y = "n_CTGs", title = NULL) +
    geom_hline(yintercept = yintercept, color = "grey50", linetype = "dashed", linewidth = 0.3) +
    scale_x_discrete(limits = levels(seu.ls.kept[[i]]$cell_type)) +
    scale_fill_manual(values = c("#B3B3B3", "#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")) +
    NoLegend() +
    ggpubr::stat_compare_means(comparisons = ls, tip.length = 0) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10))
  
  if (i %in% tail(names(seu.ls.kept), n = 1) == FALSE) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  return(p)
}
plots |> wrap_plots(nrow = 2)



## Fig. 5C-D ---------------------------------------------------------------------
# rho.df  <- foreach(i = names(ST.meta.ls), .combine = 'rbind') %do% {
#   meta <- ST.meta.ls[[i]]
#   vars <- c("e_CTGs", "n_CTGs", "n_Genes", "SingleR_malignant_score")
#   df   <- meta |> select(any_of(vars))
#   map_df(df, \(x) {
#     fit <- cor.test(x, meta$ESTIMATE_purity, method = "spearman")
#     data.frame(rho = fit$estimate, p.val = fit$p.value)
#   }) |>
#     mutate(unique_sample_id = i, var = names(df), .before = everything()) |>
#     remove_rownames()
# }
# rho.df  <- rho.df |> inner_join(ST.data, by = "unique_sample_id")
# fig.5d  <- rho.df


##### (1) C #####
# df            <- rho.df |> filter(var == "n_CTGs", with_image == "yes") |> group_by(cancer) |> slice_max(order_by = rho, n = 3) |> ungroup()
# samples       <- c("HTAN_S18-25657-Fs1U1Bp1", "HTAN_S21-46983-A3U1Bp1", "HTAN_S16-38794-A3U1Bp1", "HTAN_TB-4951U1Bp1",
#                    "SCAR_ST_000099", "Mendeley_skrx2fz79n_HCC_P3T")
# seu.ls        <- foreach(i = samples) %do% {
#   path  <- df$processed_path[df$unique_sample_id == i]
#   seu   <- readRDS(path)
#   seu$x <- log10(seu$n_CTGs + 1)
#   seu
# }
# names(seu.ls) <- samples


load("../data/fig.5c.RData")  ## c("df", "seu.ls")
min     <- sapply(seu.ls, \(x) min(x$ESTIMATE_purity)) |> min()
max     <- sapply(seu.ls, \(x) max(x$x)) |> max()
plots   <- foreach(i = names(seu.ls)) %do% {
  seu <- seu.ls[[i]]
  p1  <- SpatialFeaturePlot(seu, features = "ESTIMATE_purity", combine = FALSE)[[1]] +
    labs(fill = "Tumor Purity", subtitle = paste("Sample ID:", i), title = df$cancer[df$unique_sample_id == i]) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5, color = "#7f8c8d", size = 5)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, frame.colour = "black", ticks.colour = "black")) +
    scale_fill_gradientn(colors = c("#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"), limits = c(min, 1))
  p2  <- SpatialFeaturePlot(seu, features = "x", combine = FALSE)[[1]] +
    labs(fill = expression(log[10]("n_CTGs"))) +
    coord_fixed(ratio = 1) +
    theme(legend.position = "right", aspect.ratio = 1) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, frame.colour = "black", ticks.colour = "black")) +
    scale_fill_gradientn(colors = c("#5E4FA2", "#66C2A5", "#E6F598", "#FEE08B", "#F46D43", "#9E0142"), limits = c(0, max))
  
  p   <- p1 + p2 + plot_layout(ncol = 1) & 
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
}
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)
legends <- plot_grid(legend1, legend2, ncol = 1)


plot_grid(plotlist = plots, nrow = 1) |>
  plot_grid(legends, nrow = 1, rel_widths = c(6, 1))


##### (2) D #####
rho.df <- fig.5d
p      <- rho.df |> filter(var == "n_CTGs") |> 
  ggplot(aes(x = reorder(cancer, rho, median), y = rho)) + 
  geom_violin(fill = "#F0F0F0", color = "#404040", alpha = 0.8, linewidth = 0.5, 
              scale = "width", width = 0.7) +
  geom_jitter(aes(fill = rho), size = 2, alpha = 0.8, shape = 21, color = "grey30", stroke = 0.3,
              width = 0.1, height = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.2, color = "#2D2D2D", linewidth = 0.4) +
  scale_fill_gradient2(low = "#7F9AA8", mid = "#F0F0F0", high = "#D6604D",                        
                       midpoint = 0, limits = c(-1, 1),
                       name = expression(rho)) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.5, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#606060", linewidth = 0.4) +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = expression(paste("Spearman's ", rho)), title = "") +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.8, 0.08),
        legend.direction = "horizontal")
p



## Fig. 5E ---------------------------------------------------------------------
CNV.score.ls      <- readRDS("./Proj_CTAs/output/6. n_CTGs/ST/infercnv/CNV.score.ls.rds")
## Note: the immune score of the reference cluster should be significantly higher than others.
immune.score.diff <- foreach(i = names(CNV.score.ls), .combine = 'rbind') %do% {
  cnv.score <- CNV.score.ls[[i]]
  cnv.score$group <- ifelse(cnv.score$cnv_cluster == "normal", "ref", "obs")
  p.val     <- wilcox.test(immune_score ~ group, data = cnv.score)$p.value
  y         <- factor(cnv.score$group, levels = c("obs", "ref"))
  auc       <- pROC::auc(y, cnv.score$immune_score) |> as.numeric()
  data.frame(unique_sample_id = i, auc, p.val)
}
immune.score.diff <- immune.score.diff |> mutate(p.adj = p.adjust(p.val, method = "fdr"))
dataset.kept      <- immune.score.diff |> filter(auc > 0.6, p.adj < 0.01) |> pull(unique_sample_id)


df     <- foreach(i = names(CNV.score.ls), .combine = 'rbind') %do% {
  cnv.score            <- CNV.score.ls[[i]]
  cnv.score.by.cluster <- cnv.score |> group_by(cnv_cluster) |> summarise(cnv_score = mean(cnv_score))
  cnv.high.cluster     <- cnv.score.by.cluster |> slice_max(order_by = cnv_score, n = 2) |> pull(cnv_cluster)
  cnv.low.cluster      <- cnv.score.by.cluster |> slice_min(order_by = cnv_score, n = 2) |> pull(cnv_cluster)
  
  df1 <- cnv.score |> filter(cnv_cluster %in% cnv.high.cluster) |> summarise(cnv_score = median(cnv_score), n_CTGs = median(n_CTGs), ESTIMATE_purity = median(ESTIMATE_purity))
  df2 <- cnv.score |> filter(cnv_cluster %in% cnv.low.cluster) |> summarise(cnv_score = median(cnv_score), n_CTGs = median(n_CTGs), ESTIMATE_purity = median(ESTIMATE_purity))
  df  <- rbind(df1, df2) |> mutate(unique_sample_id = i, cnv_group = c("CNV-high", "CNV-low"))
}
fig.5e <- df |> filter(unique_sample_id %in% dataset.kept) 


p <- fig.5e |>
  ggplot(aes(x = cnv_group, y = log10(n_CTGs + 1))) + 
  geom_boxplot(aes(fill = cnv_group), width = 0.6, alpha = 0.2, linewidth = 0.3, outliers = FALSE) +
  geom_line(aes(group = unique_sample_id), color = alpha("grey75", 0.5), linewidth = 0.1, linetype = "dashed") +
  geom_point(aes(fill = cnv_group, shape = cnv_group), color = "grey30", size = 2, stroke = 0.1, alpha = 0.8,    
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  ggpubr::stat_compare_means(comparisons = list(c("CNV-high", "CNV-low")), paired = TRUE, tip.length = 0) +
  scale_fill_manual(values = c("CNV-high" = "#2171B5", "CNV-low" = "#FDB863")) +
  scale_shape_manual(values = c(21, 24)) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = expression(log[10]("median n_CTGs"))) +
  theme(legend.position = "none",  axis.text.x = element_text(angle = 30, hjust = 1), 
        panel.background = element_rect(fill = "grey98"),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3))
p

