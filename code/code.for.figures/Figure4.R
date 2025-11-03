library(tidyverse)
library(foreach)
library(doParallel)
library(ggpubr)
library(ComplexHeatmap)
library(patchwork)



load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")  ## "tumor.samples"
load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")  ## "dataset.info.10x"
load("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/CTG.score.RData")  ## c("dataset.info", "CTG.score.ls.dataset.level", "CTG.score.ls.sample.level")
source("./Proj_CTAs/code/visualization.R")


dataset.info  <- dataset.info |> 
  mutate(dataset_group = case_match(
    dataset,
    dataset.info.10x$dataset       ~ "Discovery Cohort",
    dataset.info.SmartSeq2$dataset ~ "Validation Cohort 1",
    .default = "Validation Cohort 2")) |>
  mutate(tech_group = case_match(
    tech,
    "10x"                            ~ "10x",
    "SmartSeq2"                      ~ "Smart-seq2",
    "inDrop"                         ~ "inDrop",
    "Drop-seq"                       ~ "Drop-seq",
    "Microwell array-based platform" ~ "Microwell array",
    c("Seq-Well", "Seq-Well S3")     ~ "Seq-Well"
  ))
tumor.samples <- tumor.samples |> 
  dplyr::rename(sample_group = group) |>
  mutate(dataset_sample = paste(dataset, sample, sep = "_")) |>
  left_join(dataset.info, by = "dataset")
score.colors  <- c("expression_CTGs" = "#4E79A7", "n_CTGs" = "#E15759", "n_Genes" = "#59A14F", "n_random_genes" = "#B07AA2")
tech.colors   <- c("10x" = "#1B9E77", "Smart-seq2" = "#D95F02", "inDrop" = "#7570B3", "Drop-seq" = "#E6AB02",  "Microwell array" = "#E7298A", "Seq-Well" = "#66A61E")



## Fig. 4A ---------------------------------------------------------------------
# col.df <- dataset.info |> inner_join(cancer.abbr) |> filter(tech %in% c("10x"), malignant_no == "yes", coarse_cell_type == "yes")
# col.df <- col.df |> arrange(cancer) |> column_to_rownames("dataset") |> mutate(abbr = factor(abbr, levels = cancer.abbr$abbr))
# 
# 
# median.n.CTGs <- foreach(i = rownames(col.df), .combine = 'rbind') %do% {
#   CTG.score.ls.dataset.level[[i]] |> 
#     group_by(my_coarse_cell_type) |> 
#     summarise(n = median(n_CTGs)) |> 
#     mutate(dataset             = i, 
#            my_coarse_cell_type = ifelse(my_coarse_cell_type %in% "Epithelial", "Normal Epithelial", my_coarse_cell_type))
# }
# mat           <- median.n.CTGs |> filter(my_coarse_cell_type %in% c("Other", "Skin_Melanocyte", "Unknown") == FALSE) |> pivot_wider(names_from = dataset, values_from = n) |> column_to_rownames("my_coarse_cell_type") |> as.matrix()


load("../data/fig.4a.RData")
col.ha  <- HeatmapAnnotation(Cancer = col.df$abbr, col = list(Cancer = cancer.colours), 
                            show_legend = FALSE, show_annotation_name = FALSE)
max     <- max(mat, na.rm = TRUE)
col.fun <- circlize::colorRamp2(seq(1, max), colorRampPalette(c("dodgerblue4", "peachpuff", "deeppink4"))(40))
p       <- Heatmap(mat, col = col.fun, name = "n", na_col = "white",
                   cluster_rows = TRUE, cluster_columns = FALSE,
                   show_row_names = TRUE, show_column_names = FALSE, 
                   column_split = col.df$abbr, column_title_rot = 45, border = TRUE,
                   top_annotation = col.ha)
p       <- draw(p)



## Fig. 4C ---------------------------------------------------------------------
## at least contain 10 non-malignant cells for each tumor sample (1043)
ds       <- tumor.samples |> filter(num_non_malignant >= 10) |> pull(dataset_sample) |> intersect(names(CTG.score.ls.sample.level))
auc.df   <- foreach(i = ds, .combine = 'rbind') %do% {
  vars <- c("n_CTGs", "e_CTGs", "n_Genes", "n_CTGs_control")
  df   <- CTG.score.ls.sample.level[[i]][, vars]
  y    <- factor(CTG.score.ls.sample.level[[i]]$malignant, levels = c("no", "yes"))
  auc  <- vapply(df, \(x) pROC::roc(y, x, direction = "<")$auc, 0)
  data.frame(dataset_sample = i, var = vars, auc = auc, row.names = NULL)
}
auc.df   <- auc.df |> 
  left_join(tumor.samples, by = "dataset_sample") |>
  mutate(x = case_match(
    var, 
    "e_CTGs"         ~ "expression_CTGs", 
    "n_CTGs_control" ~ "n_random_genes", 
    .default = var
  ))
auc.df$x <- factor(auc.df$x, levels = c("expression_CTGs", "n_CTGs", "n_Genes", "n_random_genes"))


df     <- auc.df |> 
  filter(dataset_group == "Discovery Cohort") |>
  left_join(cancer.abbr) |> 
  mutate(abbr = factor(abbr, levels = cancer.abbr$abbr))
fig.4c <- df


df <- fig.4c
df |> filter(var == "n_CTGs") |> group_by(abbr) |> summarise(auc = median(auc)) |> arrange(desc(auc))
df |> filter(var == "n_CTGs") |> pull(auc) |> median()
p  <- df |> filter(var == "n_CTGs") |> 
  ggplot(aes(x = reorder(abbr, auc, median), y = auc)) + 
  geom_violin(fill = "#F5F5F5", color = "#252525", alpha = 0.8, linewidth = 0.5, 
              scale = "width", width = 0.7) +
  geom_jitter(aes(fill = auc), size = 2, alpha = 0.8, shape = 21, color = "grey30", stroke = 0.3,
              width = 0.1, height = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.2, color = "#2D2D2D", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",                        
                       midpoint = 0.5, limits = c(0, 1),
                       name = "ROC-AUC") +
  scale_y_continuous(limits = c(0.15, 1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#606060", linewidth = 0.4) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "ROC-AUC", title = "Discovery Cohort") +
  theme(axis.text.x = element_text(hjust = 1, angle = 30), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.4),
        legend.position = c(0.8, 0.08),
        legend.direction = "horizontal")
p



## Fig. 4E-F -------------------------------------------------------------------
# CTG.SingleR.cnv.score.ls <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/CTG.SingleR.cnv.score.ls.rds")
# scCancer2.score.ls       <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/scCancer2.score.ls.rds")
# ikarus.score.ls          <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/ikarus.score.ls.rds")
# ikarus.samples           <- names(ikarus.score.ls)


# length(scCancer2.score.ls) < nrow(sample.info)  ## some samples were filtered when calculate CNV score
# ds              <- tumor.samples |> filter(num_non_malignant >= 10) |> pull(dataset_sample) |> intersect(names(CTG.SingleR.cnv.score.ls))
# score.ls        <- foreach(i = ds) %do% {
#   df <- inner_join(CTG.SingleR.cnv.score.ls[[i]], scCancer2.score.ls[[i]], by = "barcode")
#   if (i %in% names(ikarus.score.ls)) {
#     ikarus <- ikarus.score.ls[[i]] |> mutate(barcode = as.character(barcode))
#     df     <- inner_join(df, ikarus, by = "barcode")
#   }
#   return(df)
# }
# names(score.ls) <- ds
# vars            <- c("cnv_score", "n_CTGs", "ikarus_malignant_score", "scCancer2_malignant_score", "SingleR_malignant_score")



##### (1) E #####
# i      <- "Data_Breast_Gao2021_Breast_TNBC3"  ## a good example
# df     <- CTG.SingleR.cnv.score.ls[[i]] |> mutate(cell_type = ifelse(cell_type == "Epithelial", "Normal Epithelial", cell_type)) |> filter(cell_type != "Unknown")
# fig.4e <- df


i      <- "Data_Breast_Gao2021_Breast_TNBC3"
df     <- fig.4e
colors <- c(Malignant = "#E76F51", `Normal Epithelial` = "#2A9D8F", Stromal = "#9B59B6", Immune = "#8BB8D6", Brain_Glial = "#F4A261", "Unknown" = "#7D8B9C")
p1     <- df |>
  ggplot(aes(x = n_CTGs, y = cnv_score, color = cell_type)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 14) +
  labs(y = "CNV score", color = "Cell type") +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
p2 <- p1 %+% aes(x = SingleR_malignant_score) %+% labs(x = "SingleR-CTG malignant score")


p1 + p2 + plot_layout(ncol = 2, guides = "collect") & 
  plot_annotation(title = i, theme = theme(plot.title = element_text(hjust = 0.5))) &
  theme(legend.position = "top") &
  guides(color = guide_legend(override.aes = list(size = 2)))



##### (2) F #####
# auc.df      <- foreach(i = ds, .combine = 'rbind') %do% {
#   score.df  <- score.ls[[i]]
#   method    <- intersect(vars, names(score.df))
#   y         <- ifelse(score.df$malignant == "yes", 1, 0)
#   pr.auc    <- vapply(method, \(x) {PRROC::pr.curve(scores.class0 = score.df[, x], weights.class0 = y)$auc.integral}, 0)
#   y         <- factor(score.df$malignant, levels = c("no", "yes"))
#   roc.auc   <- vapply(method, \(x) {pROC::roc(y, score.df[, x], direction = "<")$auc}, 0)
#   
#   data.frame(dataset_sample = i, method, pr.auc, roc.auc, row.names = NULL)
# }
# auc.df      <- auc.df |> 
#   mutate(method = case_match(
#     method, 
#     "cnv_score"                 ~ "inferCNV", 
#     "n_CTGs"                    ~ "n_CTGs", 
#     "ikarus_malignant_score"    ~ "ikarus", 
#     "scCancer2_malignant_score" ~ "scCancer2", 
#     "SingleR_malignant_score"   ~ "SingleR_CTG"))
# auc.df      <- auc.df |> 
#   mutate(method = factor(method, levels = c("n_CTGs", "SingleR_CTG", "ikarus", "scCancer2", "inferCNV"))) |>
#   left_join(tumor.samples, by = "dataset_sample")



load("../data/fig.4f.RData")
create.plot <- function(df, y, title = NULL) {
  df |>
    ggplot(aes(x = method, y = .data[[y]])) +
    geom_boxplot(aes(color = method), outliers = FALSE, linewidth = 0.2) +
    geom_point(aes(color = method), size = 0.01) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
    ggpubr::stat_compare_means(comparisons = ls, tip.length = 0, paired = TRUE, group = "dataset_sample", step.increase = 0.05, size = 2.8) +
    ggpubr::stat_compare_means(data = filter(df, dataset_sample %in% ikarus.samples),
                               comparisons = list(c("SingleR_CTG", "ikarus")), tip.length = 0, paired = TRUE, group = "dataset_sample", size = 2.8) +
    theme_bw(base_size = 14) +
    labs(x = NULL, y = ifelse(grepl("roc", y), "ROC-AUC", "PR-AUC"), title = title) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = df |> filter(method == "SingleR_CTG") |> pull(get(y)) |> median(),
               linetype = "dashed", linewidth = 0.1)
}


library(aplot)
ls    <- list(c("n_CTGs", "SingleR_CTG"), c("SingleR_CTG", "scCancer2"), c("SingleR_CTG", "inferCNV"))
p1    <- create.plot(filter(auc.df, dataset_group == "Discovery Cohort"), "roc.auc", "Discovery Cohort")
p2    <- create.plot(filter(auc.df, dataset_group == "Discovery Cohort"), "pr.auc")
p3    <- create.plot(filter(auc.df, dataset_group == "Validation Cohort 1"), "roc.auc", "Validation Cohort 1") + labs(y = NULL)
p4    <- create.plot(filter(auc.df, dataset_group == "Validation Cohort 1"), "pr.auc") + labs(y = NULL)
p5    <- create.plot(filter(auc.df, dataset_group == "Validation Cohort 2"), "roc.auc", "Validation Cohort 2") + labs(y = NULL)
p6    <- create.plot(filter(auc.df, dataset_group == "Validation Cohort 2"), "pr.auc") + labs(y = NULL)


(p1 |> insert_right(p3) |> insert_right(p5) |> as.patchwork()) /
  (p2 |> insert_right(p4) |> insert_right(p6) |> as.patchwork())



## Fig. 4G----------------------------------------------------------------------
# load("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/running.time.input.RData")
# time.df <- readRDS("./Proj_CTAs/output/6. n_CTGs/scRNA-seq/time.df.rds")
# time.df <- time.df |> 
#   left_join(args.df, by = "dataset_sample") |> 
#   select(dataset_sample, method, time, n_cell, n_cell_group) |>
#   mutate(n_cell_group = factor(paste0("~", n_cell_group), levels = c("~1k", "~3k", "~5k", "~8k", "~10k")))
# fig.4g  <- time.df


time.df <- fig.4g
colors  <- c("SingleR" = "#A3A500", "infercnv" = "#E76BF3")
m       <- "SingleR"
df      <- filter(time.df, method == m) |> mutate(time = as.numeric(time))
df.summ <- df |> group_by(n_cell_group) |> summarise(mean = mean(time))
p1      <- ggplot(df, aes(x = n_cell_group, y = time)) +
  geom_line(data = df.summ, aes(x = n_cell_group, y = mean, group = 1),
            color = colors[m], linewidth = 1) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 2, fill = colors[m], shape = 21, color = "black") +
  geom_point(data = df.summ, aes(x = n_cell_group, y = mean),
             shape = 21, size = 4, fill = colors[m], color = "white", stroke = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.15, color = colors[m], linewidth = 0.8, alpha = 0.8) +
  labs(x = "Number of cells", y = "Runtime (seconds)", title = "SingleR-CTG") +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
        plot.title = element_text(hjust = 0.5)) 


m       <- "infercnv"
df      <- filter(time.df, method == m) |> mutate(time = as.numeric(time) / 60)
df.summ <- df |> group_by(n_cell_group) |> summarise(mean = mean(time))
p2      <- ggplot(df, aes(x = n_cell_group, y = time)) +
  geom_line(data = df.summ, aes(x = n_cell_group, y = mean, group = 1),
            color = colors[m], linewidth = 1) +
  geom_jitter(width = 0.1, alpha = 0.25, size = 2, fill = colors[m], shape = 21, color = "black") +
  geom_point(data = df.summ, aes(x = n_cell_group, y = mean),
             shape = 21, size = 4, fill = colors[m], color = "white", stroke = 0.8) +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.15, color = colors[m], linewidth = 0.8, alpha = 0.8) +
  labs(x = "Number of cells", y = "Runtime (minutes)", title = m) +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
        plot.title = element_text(hjust = 0.5)) 


p1 / p2

