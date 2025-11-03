library(tidyverse)
library(foreach)
library(magrittr)
library(ggpattern)
library(patchwork)
library(Seurat)



## Fig. 7A ---------------------------------------------------------------------
# load("./Proj_CTAs/output/8. clinical relevance/survival.RData")
# df       <- foreach(i = names(HR.df.ls), .combine = 'rbind') %do% {
#   HR.df.ls[[i]] |> 
#     filter(!is.na(HR), CTG_group %in% c("CTG in this cancer", "not CTG")) |>
#     mutate(p.adj         = p.adjust(p.val, method = "fdr"),
#            gene_group    = case_when(
#              HR > 1 & p.adj < 0.05 ~ "HR > 1, p.adj < 0.05",
#              HR < 1 & p.adj < 0.05 ~ "HR < 1, p.adj < 0.05",
#              .default = "not correlated"),
#            is_correlated = ifelse(gene_group %in% c("HR > 1, p.adj < 0.05", "HR < 1, p.adj < 0.05"), "yes", "no"),
#            CTG_group     = ifelse(CTG_group == "not CTG", "Non_CTG", "CTG") |> factor(levels = c("Non_CTG", "CTG")),
#            cancer_type   = i)
# }
# df       <- df |> mutate(cancer_type = factor(cancer_type, levels = names(HR.df.ls)))
# fig.7a   <- df


df       <- fig.7a
format_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1), round(p, 3))
}
pval     <- df |> 
  group_by(cancer_type) |>
  group_modify(~ {
    tb    <- .x %$% table(is_correlated, CTG_group)
    p.val <- if (identical(dim(tb), as.integer(c(2, 2)))) {
      fisher.test(tb)$p.value
    } else {
      NA
    }
    data.frame(p.val) |> mutate(label = paste("p =", format_p(p.val)))
  }) |>
  mutate(label = ifelse(label == "p = NA", "", label))


p <- ggplot(df, aes(x = CTG_group, fill = gene_group)) +
  geom_bar(position = "fill", width = 0.6, color = "#F5F5F5", linewidth = 0.02) +
  facet_wrap(~ cancer_type, ncol = 13, scales = "free_x") +
  geom_text(data = pval, aes(x = 1.5, y = 1.05, label = label), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(
    name   = NULL,
    values = c("HR > 1, p.adj < 0.05" = "#C23522", "HR < 1, p.adj < 0.05" = "#2B5876", "not correlated" =  "grey90"),  ## "#F5F5F5"
    breaks = c("HR > 1, p.adj < 0.05", "HR < 1, p.adj < 0.05")) +
  scale_y_continuous(
    name = "Proportion of survival-associated genes",
    labels = scales::percent_format(),
    breaks = seq(0, 1, length.out = 6),
    expand = expansion(mult = c(0, 0.05))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "bottom")
p



## Fig. 7B ---------------------------------------------------------------------
# source("./Proj_CTAs/code/visualization.R")
# load("./Proj_CTAs/output/1. data wrangling/datasets.summary.RData")
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# load("./Proj_CTAs/output/1. data wrangling/pure.tumor.expr.RData")
# 
# 
# LM22.percent <- readRDS("./Proj_CTAs/output/2. cell annotation/LM22.percent.rds")
# samples.ls   <- LM22.percent |>
#   distinct(cancer, dataset, tech, dataset_sample, sample_group) |>
#   filter(dataset %in% dataset.info.10x$dataset,
#          dataset_sample %in% with(tumor.samples, paste(dataset, sample, sep = "_")),
#          sample_group == "primary") %$%
#   split(dataset_sample, cancer)
# percent      <- LM22.percent |> dplyr::rename(immune_cell = LM22)
# cells        <- levels(LM22.percent$LM22)
# 
# 
# hallmarks        <- msigdbr::msigdbr(category = "H") %$% split(gene_symbol, gs_name)
# names(hallmarks) <- names(hallmarks) |> str_sub(start = 10)
# immune.hallmarks <- c("ALLOGRAFT_REJECTION", "COMPLEMENT", "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE",
#                       "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE")
# 
# 
# ssgsea.ls        <- foreach(i = cancers) %do% {
#   s      <- samples.ls[[i]]
#   expr   <- pure.tumors.10x.cpm[[i]]
#   ssgsea <- GSVA::gsva(expr          = expr[, s], 
#                        gset.idx.list = c(list(CTGs = pan.sample.CTGs.ls[[i]]), hallmarks),
#                        method        = "ssgsea",
#                        abs.ranking   = FALSE,
#                        ssgsea.norm   = TRUE,
#                        BPPARAM       = BiocParallel::MulticoreParam(40))
# }
# names(ssgsea.ls) <- cancers 
# 
# 
# ssgsea.df         <- foreach(i = cancers, .combine = 'rbind') %do% {
#   ssgsea.ls[[i]] |> t() |> as.data.frame() |> rownames_to_column("dataset_sample") |> mutate(cancer = i)
# }
# ssgsea.percent.df <- ssgsea.df |> 
#   left_join(percent, by = c("cancer", "dataset_sample")) |> 
#   left_join(cancer.abbr) |> 
#   mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) |>
#   filter(n_immune >= 100)
# fig.7b            <- ssgsea.percent.df



ssgsea.percent.df <- fig.7b
n.pri.samples     <- ssgsea.percent.df |> 
  distinct(abbr, dataset_sample) |> 
  dplyr::count(abbr) |>
  mutate(x = paste0(abbr, "\n(n = ", n, ")"),
         x = factor(x, levels = x))
cells             <- levels(ssgsea.percent.df$immune_cell)


cell.rho.df <- ssgsea.percent.df |>
  group_by(abbr) |>
  group_modify(~ {
    foreach(i = cells, .combine = 'rbind') %do% {
      fit        <- .x |> filter(immune_cell == i) %$% cor.test(CTGs, percentage_immune, method = "spearman")
      percentage <- .x |> filter(immune_cell == i) |> pull(percentage_immune) |> median()
      data.frame(immune_cell = i, percentage, rho = fit$estimate, p.val = fit$p.value, row.names = NULL)
    }
  }) |>
  mutate(immune_cell = factor(immune_cell, levels = cells))
cells.kept  <- c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated", "T cells follicular helper", "T cells regulatory (Tregs)", 
                 "NK cells resting", "NK cells activated", 
                 "Macrophages M0", "Macrophages M1", "Macrophages M2")
rho.df      <- cell.rho.df |> left_join(n.pri.samples) |> filter(n >= 10, immune_cell %in% cells.kept)
na.df       <- rho.df |> filter(is.na(rho))
valid.df    <- rho.df |> filter(!is.na(rho))
p.sig       <- 0.05
text.df     <- valid.df |>
  filter(p.val <= p.sig) |>
  mutate(label      = sprintf("%.2f", rho),
         text_color = ifelse(rho > 0, "white", "black"),
         sig_star   = case_when(
           p.val <= 0.0001 ~ "****",
           p.val <= 0.001 ~ "***",
           p.val <= 0.01 ~ "**",
           p.val <= 0.05 ~ "*")
  )
p1          <- ggplot() +
  geom_tile(data = valid.df, aes(x = x, y = immune_cell, fill = rho),
            color = "white", linewidth = 0.2) +
  geom_tile_pattern(data = na.df, aes(x = x, y = immune_cell),
                    color = "white", fill = "grey90",
                    pattern = "crosshatch", pattern_color = "grey60", pattern_fill = "grey90",
                    pattern_alpha = 0.5, pattern_density = 0.25, pattern_spacing = 0.015) +
  ## add black box
  geom_tile(data = filter(valid.df, p.val <= p.sig), aes(x = x, y = immune_cell),
            color = "black", fill = NA, linewidth = 0.25, width = 0.95, height = 0.95) +
  ## add rho
  geom_text(data = text.df, aes(x = x, y = immune_cell, label = label, color = text_color),
            size = 3.0, show.legend = FALSE) +
  ## add pval
  geom_text(data = text.df, aes(x = x, y = immune_cell, label = sig_star),
            color = "black", size = 2.5, vjust = -0.25, hjust = 0.5, show.legend = FALSE) +
  scale_fill_gradient2(low = "#7F9AA8", mid = "#F0F0F0", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = expression(rho),
                       na.value = "transparent",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_identity() +  ## directly using the colors recorded in "text_color"
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 8), panel.grid = element_blank())


pathway.rho.df <- ssgsea.percent.df |>
  group_by(abbr) |>
  group_modify(~ {
    foreach(i = names(ssgsea.percent.df)[3:52], .combine = 'rbind') %do% {
      fit <- cor.test(.x$CTGs, .x[[i]], method = "spearman")
      data.frame(pathway = i, rho = fit$estimate, p.val = fit$p.value, row.names = NULL)
    }
  })
pathways.kept  <- c("INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE")
rho.df         <- pathway.rho.df |> left_join(n.pri.samples) |> 
  filter(n >= 10, pathway %in% pathways.kept) |>
  mutate(pathway = pathway |> str_replace_all("_", " ") |> str_to_sentence())
na.df          <- rho.df |> filter(is.na(rho))
valid.df       <- rho.df |> filter(!is.na(rho))
p.sig          <- 0.05
text.df        <- valid.df |>
  filter(p.val <= p.sig) |>
  mutate(label      = sprintf("%.2f", rho),
         text_color = ifelse(rho > 0, "white", "black"),
         sig_star   = case_when(
           p.val <= 0.0001 ~ "****",
           p.val <= 0.001 ~ "***",
           p.val <= 0.01 ~ "**",
           p.val <= 0.05 ~ "*")
  )
p2             <- ggplot() +
  geom_tile(data = valid.df, aes(x = x, y = pathway, fill = rho),
            color = "white", linewidth = 0.2) +
  geom_tile_pattern(data = na.df, aes(x = x, y = pathway),
                    color = "white", fill = "grey90",
                    pattern = "crosshatch", pattern_color = "grey60", pattern_fill = "grey90",
                    pattern_alpha = 0.5, pattern_density = 0.25, pattern_spacing = 0.015) +
  ## add black box
  geom_tile(data = filter(valid.df, p.val <= p.sig), aes(x = x, y = pathway),
            color = "black", fill = NA, linewidth = 0.25, width = 0.95, height = 0.95) +
  ## add rho
  geom_text(data = text.df, aes(x = x, y = pathway, label = label, color = text_color),
            size = 3.0, show.legend = FALSE) +
  ## add pval
  geom_text(data = text.df, aes(x = x, y = pathway, label = sig_star),
            color = "black", size = 2.5, vjust = -0.25, hjust = 0.5, show.legend = FALSE) +
  scale_fill_gradient2(low = "#7F9AA8", mid = "#F0F0F0", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1), name = expression(rho),
                       na.value = "transparent",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_color_identity() +  ## directly using the colors recorded in "text_color"
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())


p1 + p2 +
  plot_layout(ncol = 1, guides = "collect", heights = c(11, 2)) &
  theme(legend.position = "bottom")



## Fig. 7C ---------------------------------------------------------------------
library(survival)
library(survminer)


# source("./Proj_CTAs/code/visualization.R")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")
# CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# 
# 
# tcga.survival <- data.table::fread("./Datas/UCSC_Xena/PANCAN/Survival_SupplementalTable_S1_20171025_xena_sp")
# tcga.survival <- tcga.survival |> 
#   rename(patient = `_PATIENT`, cancer_abbr = `cancer type abbreviation`) |> 
#   mutate(sample_type_code = as.numeric(str_sub(sample, start = 14))) |>
#   filter(sample_type_code < 10) |>
#   slice_min(order_by = sample_type_code, n = 1, by = patient)
# 
# 
# target.survival        <- data.table::fread("./Datas/UCSC_Xena/TARGET/TARGET_donor_allprojects_transfer_to_sample.gz")
# names(target.survival) <- gsub("^_", "", names(target.survival))
# stopifnot(identical(target.survival$OS_IND, target.survival$EVENT))
# stopifnot(identical(target.survival$OS, target.survival$TIME_TO_EVENT))
# target.survival$xena_sample |> duplicated() |> table()  ## no duplicated samples
# 
# 
# survival.df <- rbind(
#   tcga.survival |> select(sample, OS, OS.time),
#   target.survival |> mutate(OS = ifelse(EVENT == "Dead", 1L, 0L)) |> select(sample = xena_sample, OS, OS.time = TIME_TO_EVENT)
# )
# survival.df <- survival.df |> filter(!is.na(OS), !is.na(OS.time))
# 
# 
# ## survival + expression data
# hgnc          <- data.table::fread("./Datas/HGNC/non_alt_loci_set.txt")
# hgnc.protein  <- data.table::fread("./Datas/HGNC/gene_with_protein_product.txt")
# protein.genes <- gene.info |> filter(ensembl %in% hgnc.protein$ensembl_gene_id, id %in% rownames(log2.tpm)) |> pull(id)
# 
# 
# merge.df.ls        <- foreach(i = cancer.abbr$cancer) %do% {
#   CTGs.ensembl <- CTGs.df |> filter(symbol %in% pan.sample.CTGs.ls[[i]]) |> pull(ensembl)
#   s            <- tumor.sample.ls[[i]]$sample
#   ssgsea       <- GSVA::gsva(expr          = log2.tpm[protein.genes, s], 
#                              gset.idx.list = list(CTGs = gene.info |> filter(ensembl %in% CTGs.ensembl) |> pull(id)),
#                              method        = "ssgsea",
#                              abs.ranking   = FALSE,
#                              ssgsea.norm   = TRUE,
#                              BPPARAM       = BiocParallel::MulticoreParam(40))
#   expr.df      <- t(ssgsea) |> as.data.frame() |> rownames_to_column("sample")
#   
#   merge.df     <- survival.df |> 
#     inner_join(expr.df) |> 
#     mutate(group = ifelse(CTGs > median(CTGs), "High", "Low"))
# }
# names(merge.df.ls) <- cancer.abbr$abbr
# fig.7c             <- merge.df.ls


merge.df.ls  <- fig.7c
i            <- "HCC"
df           <- merge.df.ls[[i]]
surv.obj     <- Surv(time = df$OS.time/30, event = df$OS)
KM.curve     <- survfit(surv.obj ~ group, data = df) 
log.rank.fit <- survdiff(surv.obj ~ group, data = df)  
coex.fit     <- coxph(surv.obj ~ CTGs, data = df)


format_p <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1), round(p, 3))
}
log.rank.p <- log.rank.fit$pvalue |> format_p()
cox.summa  <- summary(coex.fit)
HR         <- cox.summa$coefficients[, "exp(coef)"] |> round(digits = 1)
HR.ci      <- cox.summa$conf.int[, c("lower .95", "upper .95")] |> round(digits = 1)
HR.p       <- cox.summa$logtest["pvalue"] |> format_p()


p <- ggsurvplot(
  fit = KM.curve,
  data = df,
  pval = FALSE, 
  risk.table = FALSE,
  conf.int = TRUE,
  palette = c("#E64B35", "#357EBD"),
  legend.title = "CTGs Expression",
  legend.labs = c(paste0("High (n = ", sum(df$group == "High"), ")"), paste0("Low (n = ", sum(df$group == "Low"), ")")),
  title = i,
  xlab = "Time (months)",
  ylab = "Overall Survival",
  break.time.by = 50) 
p <- p$plot + 
  annotate("text", x = c(0, 0, 0), y = c(0.25, 0.20, 0.15), hjust = 0,
           label = c(paste("Log-rank P =", log.rank.p),
                     paste0("HR = ", HR, " (", HR.ci[1], " - ", HR.ci[2], ")"),
                     paste("HR P =", HR.p))) +
  theme_classic2(base_size = 12) + 
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.1),
    plot.title = element_text(hjust = 0.5))
p



## Fig. 7D ---------------------------------------------------------------------
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")
# CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# 
# 
# hgnc          <- data.table::fread("./Datas/HGNC/non_alt_loci_set.txt")
# hgnc.protein  <- data.table::fread("./Datas/HGNC/gene_with_protein_product.txt")
# protein.genes <- gene.info |> filter(ensembl %in% hgnc.protein$ensembl_gene_id, id %in% rownames(log2.tpm)) |> pull(id)
# 
# 
# i            <- "Liver Cancer (HCC)"
# CTGs.ensembl <- CTGs.df |> filter(symbol %in% pan.sample.CTGs.ls[[i]]) |> pull(ensembl)
# s            <- tumor.sample.ls[[i]]$sample
# ssgsea       <- GSVA::gsva(expr          = log2.tpm[protein.genes, s], 
#                            gset.idx.list = list(CTGs = gene.info |> filter(ensembl %in% CTGs.ensembl) |> pull(id)),
#                            method        = "ssgsea",
#                            abs.ranking   = FALSE,
#                            ssgsea.norm   = TRUE,
#                            BPPARAM       = BiocParallel::MulticoreParam(40))
# expr.df      <- t(ssgsea) |> as.data.frame() |> 
#   rownames_to_column("sample") |>
#   mutate(group = ifelse(CTGs > median(CTGs), "CTG-high", "CTG-low"))
# 
# 
# LIHC.sample.info <- data.table::fread("./Datas/UCSC_Xena/PANCAN/TCGA.LIHC.sampleMap_LIHC_clinicalMatrix")
# LIHC.sample.info <- LIHC.sample.info |>
#   select(sample = sampleID, pathologic_stage) |>
#   mutate(stage = case_match(
#     pathologic_stage,
#     "Stage I"                                                ~ "Stage I",
#     "Stage II"                                               ~ "Stage II",
#     c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
#     c("Stage IV", "Stage IVA", "Stage IVB")                  ~ "Stage IV"
#   )) |>
#   filter(!is.na(stage))
# fig.7d           <- inner_join(expr.df, LIHC.sample.info, by = "sample")


ls <- lapply(c("Stage II", "Stage III", "Stage IV"), c, "Stage I")
p  <- fig.7d |> 
  ggplot(aes(x = stage, y = CTGs)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.1, color = "#2E5984", size = 1, alpha = 0.6) +
  ggpubr::stat_compare_means(comparisons = ls, tip.length = 0) +
  theme_classic(base_size = 12) +
  labs(x = "Pathologic stage", y = "CTG expression", title = "TCGA-LIHC") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), plot.title = element_text(hjust = 0.5))
p



## Fig. 7E ---------------------------------------------------------------------
# log.tpm <- data.table::fread("./Proj_CTAs/input/scRNA-seq datasets/HCC/Data_Liver_Sun2021_log2TPM/HCC_log_tpm_expression_matrix.txt.gz")
# log.tpm <- log.tpm |> column_to_rownames("gene") |> as.matrix() |> as("dgCMatrix")
# meta    <- data.table::fread("./Proj_CTAs/input/scRNA-seq datasets/HCC/Data_Liver_Sun2021_log2TPM/HCC_cell_metadata.txt") %>% .[-1, ]
# stopifnot(all(meta$name == meta$sample_name))
# stopifnot(all(meta$name == colnames(log.tpm)))
# meta    <- meta |> select(-sample_name) |> mutate(sample = str_sub(name, end = -6)) |> column_to_rownames("name")
# 
# 
# seu                 <- CreateSeuratObject(counts = log.tpm, meta.data = meta)
# seu[["RNA"]]@counts <- matrix(ncol = 0, nrow = 0)
# 
# 
# immune.ref               <- readRDS("./Proj_CTAs/input/SingleR/monaco_immune.rds")
# blueprint.ref            <- readRDS("./Proj_CTAs/input/SingleR/blueprint_encode.rds")
# predictions              <- SingleR::SingleR(test = seu[["RNA"]]@data, ref = blueprint.ref, labels = blueprint.ref$label.main)
# seu$blueprint_label_main <- predictions$pruned.labels
# predictions              <- SingleR::SingleR(test = seu[["RNA"]]@data, ref = immune.ref, labels = immune.ref$label.main)
# seu$immune_label_main    <- predictions$pruned.labels
# predictions              <- SingleR::SingleR(test = seu[["RNA"]]@data, ref = immune.ref, labels = immune.ref$label.fine)
# seu$immune_label_fine    <- predictions$pruned.labels
# ## check the annotation the author provided
# pheatmap::pheatmap(table(seu$cell_type, seu$blueprint_label_main), scale = "row")
# 
# 
# Immune <- c("C0_Tcell", "C1_Tcell", "C3_Tcell", "C5_Tcell", "C19_Tcell", 
#             "C4_NK", "C7_NK", "C6_Bcell", "C21_pDC", "C22_Plasma", 
#             "C2_Mye.", "C8_Mye.", "C11_Mye.", "C15_Mey.", "C20_Mye")
# Tumor  <- c("C9_Tumor", "C10_Tumor", "C12_Tumor", "C13_Tumor", "C14_Tumor", "C16_Tumor")
# Other  <- c("C17_Endo.", "C18_Epi.", "C23_HSC")
# Tcell  <- c("C0_Tcell", "C1_Tcell", "C3_Tcell", "C5_Tcell", "C19_Tcell")
# 
# 
# seu$malignant <- ifelse(seu$cell_type %in% Tumor, "yes", "no")
# seu$Treg      <- ifelse(seu$cell_type %in% Tcell & seu$immune_label_fine %in% "T regulatory cells", "yes", "no")
# VlnPlot(seu, features = c("CD4", "IL2RA", "FOXP3"), group.by = "Treg")
# 
# 
# Treg.percent      <- seu[[]] |> group_by(sample) |> filter(cell_type %in% Immune) |> summarise(percent = mean(Treg == "yes"))
# seu.mal           <- subset(seu, malignant == "yes" & tissue_source == "Tumor")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# HCC.CTGs          <- pan.sample.CTGs.ls$`Liver Cancer (HCC)`
# seu.mal           <- AddModuleScore(seu.mal, list(HCC.CTGs))
# seu.mal$CTG_score <- seu.mal$Cluster1
# seu.mal$Cluster1  <- NULL
# fig.7e            <- seu.mal[[]] |> group_by(sample) |> summarise(x = median(CTG_score)) |> left_join(Treg.percent) |> left_join(distinct(seu.mal[[]], sample, HCC_type))


p  <- fig.7e |>
  mutate(group = ifelse(x > median(x), "CTG-high", "CTG-low")) %>%
  ggplot(aes(x = group, y = percent)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6, color = "#E63946") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, color = "#2B5F87", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "errorbar", aes(ymin = after_stat(y), ymax = after_stat(y)),
               width = 0.3, color = "#2B5F87", linewidth = 0.8) +
  ggpubr::stat_compare_means(comparisons = list(c("CTG-high", "CTG-low")), tip.length = 0) +
  labs(x = NULL, y = "Treg Cell Abundance (% of Immune Cells)", title = "Sun et al.") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
p



## Fig. 7F ---------------------------------------------------------------------
genes <- c("GAGE1", "PAGE1", "PAGE5", "GTSF1",    ## CTGs
           "BIRC5", "TOP2A", "UBE2C",             ## cell cycle
           "EPCAM", "MDK", "HMGA1",               ## stemness
           "H2AFZ", "HIST1H4C", "HMGN2", "PTMA",  ## chromatin structure and transcriptional regulation; HMGN2 (Chromatin Organization)
           "NQO1", "FABP5",
           
           "CD74", "HLA-DRA", "C3",               ## immune
           "INSIG1",                              ## cholesterol synthesis
           "SULT2A1",
           "CYP27A1", "CYP2A6", "CYP2C9", "CYP2D6", "CYP3A4", "CYP3A5",
           "APOB", "APOA5", "APOC1", "APOC2", "APOC3")
genes <- c("GAGE1", "PAGE1", "PAGE5", "GTSF1",    ## CTGs
           "BIRC5", "TOP2A",                      ## cell cycle
           "EPCAM", "MDK", "HMGA1",               ## stemness
           "H2AFZ", "HIST1H4C", "HMGN2",          ## chromatin structure and transcriptional regulation; HMGN2 (Chromatin Organization)
           "NQO1", "FABP5",
           
           "CD74", "HLA-DRA", "C3",               ## immune
           "INSIG1",                              ## cholesterol synthesis
           "SULT2A1",
           "CYP27A1", "CYP2A6", "CYP2D6",
           "APOB", "APOC1", "APOC2", "APOC3")


df    <- fig.7f
fit   <- cor.test(df$x, df$y, method = "spearman")
label <- paste0("R = ", round(fit$estimate, digits = 2), ", p < ", "2.2e-16")
p     <- ggplot(df, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3, linetype = "dotted") +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3, linetype = "dotted") +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linewidth = 0.4, linetype = "dashed") +
  geom_point(aes(color = group, alpha = -log10(p_adj_max)), 
             size = 2, shape = 16) +
  scale_color_manual(values = c("consistently up-regulated" = "#E64B35", "consistently down-regulated" = "#3182BD", "others" = "grey70")) +
  scale_alpha_continuous(name = expression(-log[10](p.adj[max]))) +
  geom_point(data = filter(df, gene %in% genes),
             color = "firebrick", size = 3, shape = 21) +
  labs(x     = expression("log"[2]*"FC (Lu et al.)"),
       y     = expression("log"[2]*"FC (Sun et al.)"),
       color = "Group") +
  theme_classic(base_size = 10) +
  theme(panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line       = element_blank(),
        axis.ticks      = element_line(color = "black", linewidth = 0.3),
        axis.text       = element_text(color = "black", size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         alpha = guide_legend(override.aes = list(color = "grey30"))) +
  geom_rect(xmin = 1, xmax = Inf, ymin = 1, ymax = Inf, 
            fill = NA, color = "#E64B35", linewidth = 0.15, linetype = "dashed") +
  geom_rect(xmin = -Inf, xmax = -1, ymin = -Inf, ymax = -1, 
            fill = NA, color = "#3182BD", linewidth = 0.15, linetype = "dashed") +
  ggrepel::geom_text_repel(data = filter(df, gene %in% genes),
                           aes(label = gene), color = "black", segment.color = "grey30", size = 2.8,
                           box.padding = 0.4, min.segment.length = 0.1, max.overlaps = 50, 
                           force = 0.5, nudge_x = 0.15, nudge_y = 0.1, fontface = "italic") +
  annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.1, vjust = 2) +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1,
        legend.position = c(0.15, 0.78))
p



## Fig. 7G ---------------------------------------------------------------------
load("../data/fig.7g.RData")
i  <- "CytoTRACE2_Score"
ls <- list(c("Low", "Middle"), c("Low", "High"))


seu.mal    <- seu.Lu
y.max      <- max(seu.mal[[]][[i]])
yintercept <- seu.mal[[]] |> group_by(group) |> summarise(x = median(get(i))) |> pull(x) |> min()
p1         <- VlnPlot(seu.mal, features = i, pt.size = 0, y.max = y.max * (1 + length(ls) * 0.15)) + 
  stat_summary(fun = median, geom = "point", size = 4, colour = "white", shape = 18) +
  scale_fill_manual(values = c("#D4B9DA", "#C994C7", "#DF65B0")) +
  labs(x = NULL, y = gsub("_", " ", i), title = NULL) +
  geom_hline(yintercept = yintercept, color = "grey30", linetype = "dashed") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = ls, tip.length = 0)


seu.mal    <- seu.Sun
y.max      <- max(seu.mal[[]][[i]])
yintercept <- seu.mal[[]] |> group_by(group) |> summarise(x = median(get(i))) |> pull(x) |> min()
p2         <- VlnPlot(seu.mal, features = i, pt.size = 0, y.max = y.max * (1 + length(ls) * 0.15)) + 
  stat_summary(fun = median, geom = "point", size = 4, colour = "white", shape = 18) +
  scale_fill_manual(values = c("#D4B9DA", "#C994C7", "#DF65B0")) +
  labs(x = NULL, y = gsub("_", " ", i), title = NULL) +
  geom_hline(yintercept = yintercept, color = "grey30", linetype = "dashed") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  ggpubr::stat_compare_means(comparisons = ls, tip.length = 0)


p1 / p2



## Fig. 7H ---------------------------------------------------------------------
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")
# CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")
# 
# 
# hgnc          <- data.table::fread("./Datas/HGNC/non_alt_loci_set.txt")
# hgnc.protein  <- data.table::fread("./Datas/HGNC/gene_with_protein_product.txt")
# protein.genes <- gene.info |> filter(ensembl %in% hgnc.protein$ensembl_gene_id, id %in% rownames(log2.tpm)) |> pull(id)
# 
# 
# i            <- "Liver Cancer (HCC)"
# CTGs.ensembl <- CTGs.df |> filter(symbol %in% pan.sample.CTGs.ls[[i]]) |> pull(ensembl)
# s            <- tumor.sample.ls[[i]]$sample
# ssgsea       <- GSVA::gsva(expr          = log2.tpm[protein.genes, s], 
#                            gset.idx.list = list(CTGs = gene.info |> filter(ensembl %in% CTGs.ensembl) |> pull(id)),
#                            method        = "ssgsea",
#                            abs.ranking   = FALSE,
#                            ssgsea.norm   = TRUE,
#                            BPPARAM       = BiocParallel::MulticoreParam(40))
# expr.df      <- t(ssgsea) |> as.data.frame() |> 
#   rownames_to_column("sample") |>
#   mutate(group = ifelse(CTGs > median(CTGs), "CTG-high", "CTG-low"))
# 
# 
# stem.RNAexp <- read.table("./Datas/UCSC_Xena/PANCAN/StemnessScores_RNAexp_20170127.2.tsv.gz", header = TRUE, row.names = 1)
# stem.RNAexp <- stem.RNAexp |> t() |> as.data.frame() |> rownames_to_column("sample") |> mutate(sample = gsub("\\.", "-", sample)) 
# fig.7h      <- inner_join(expr.df, stem.RNAexp, by = "sample")


p <- fig.7h |>
  ggpubr::ggscatter(x = "CTGs", y = "RNAss", color = "#2E5984", size = 1, alpha = 0.6,
                    add = "reg.line", add.params = list(color = "#C85200", fill = "gray"),
                    conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "spearman", size = 3.5)) + 
  labs(x = "CTG expression", y = "Stemness score (mRNAsi)", title = "TCGA-LIHC") +
  theme_classic(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1)
p
