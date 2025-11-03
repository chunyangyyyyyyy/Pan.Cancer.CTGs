library(tidyverse)
library(patchwork)


source("./Proj_CTAs/code/visualization.R")



## Fig. 1C ---------------------------------------------------------------------
# library(foreach)
# library(doParallel)
# 
# 
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# 
# 
# cl               <- makeCluster(20)
# registerDoParallel(cl)
# ## genes for each dataset
# genes.ls         <- foreach(dataset = dataset.info.10x$dataset, .packages = 'Seurat') %dopar% {
#   prefix <- "./Proj_CTAs/input/scRNA-seq datasets/normalized Seurat object/"
#   file   <- paste0(dataset, ".rds")
#   seu    <- readRDS(paste0(prefix, file))
#   rownames(seu)
# }
# names(genes.ls)  <- dataset.info.10x$dataset
# 
# 
# n.sample <- dataset.info.10x |> group_by(cancer) |> summarise(n_tumor = sum(n_tumor))
# n.CTG    <- data.frame(n_CTG = sapply(pan.sample.CTGs.ls, length)) |> rownames_to_column("cancer")
# n.genes  <- foreach(i = df$cancer, .combine = 'rbind') %do% {
#   datasets <- dataset.info.10x |> filter(cancer == i) |> pull(dataset)
#   data.frame(cancer = i, n_gene = length(Reduce(intersect, genes.ls[datasets])))
# }
# df       <- Reduce(inner_join, list(n.sample, n.CTG, n.genes))
# df       <- df |> left_join(cancer.abbr) |> mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) 
# fig.1c   <- df


df    <- fig.1c
upper <- quantile(df$n_CTG, 0.75) + 1.5 * IQR(df$n_CTG)
lower <- quantile(df$n_CTG, 0.25) - 1.5 * IQR(df$n_CTG)
p1    <- ggplot(df, aes(x = n_tumor, y = n_CTG)) + 
  geom_point(aes(color = abbr), size = 2.5) +
  scale_color_manual(values = cancer.colours) +
  geom_point(data = filter(df, abbr %in% c("ccRCC")), color = "#D53E4F", fill = NA, size = 3.5, shape = 21, stroke = 0.8) +
  ggrepel::geom_text_repel(aes(label = abbr), color = "#252525", size = 3.8, box.padding = 0.3, segment.color = "#737373") +
  geom_hline(yintercept = c(lower, upper), linetype = "dashed", color = "#E7298A", linewidth = 0.3, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "Number of tumor samples", y = "Number of identified CTGs") +
  theme(aspect.ratio = 1, legend.position = "none")
p2    <- ggplot(df, aes(x = n_gene, y = n_CTG)) + 
  geom_point(aes(color = abbr), size = 2.5) + 
  scale_color_manual(values = cancer.colours) +
  geom_point(data = filter(df, abbr %in% c("ccRCC")), color = "#D53E4F", fill = NA, size = 3.5, shape = 21, stroke = 0.8) +
  ggrepel::geom_text_repel(aes(label = abbr), color = "#252525", size = 3.8, box.padding = 0.3, segment.color = "#737373") +
  geom_hline(yintercept = c(lower, upper), linetype = "dashed", color = "#E7298A", linewidth = 0.3, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(x = "Number of total genes", y = "Number of identified CTGs") +
  theme(aspect.ratio = 1, legend.position = "none")
p1 | p2



## Fig. 1D ---------------------------------------------------------------------
# library(foreach)
# 
# 
# load("./Proj_CTAs/output/1. data wrangling/selected.cancers.RData")
# load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
# load("./Proj_CTAs/output/3. define CTGs/pct.summary.RData")
# load("./Proj_CTAs/output/1. data wrangling/pure.tumor.expr.RData")
# 
# 
# CTGs.ls        <- foreach(i = names(pct.summary.each.cancer.SmartSeq2)) %do% {
#   intersect(pan.sample.CTGs.ls[[i]], pct.summary.each.cancer.SmartSeq2[[i]]$CTG)
# }
# names(CTGs.ls) <- names(pct.summary.each.cancer.SmartSeq2)
# CTGs.tpm       <- foreach(i = names(CTGs.ls), .combine = 'rbind') %do% {
#   pure.tumors.SmartSeq2.tpm[[i]][CTGs.ls[[i]], ] |>
#     as.data.frame() |>
#     rownames_to_column("CTG") |>
#     pivot_longer(-CTG) |>
#     group_by(CTG) |>
#     summarise(median_tpm = median(value)) |>
#     ungroup() |>
#     mutate(cancer = i)
# }
# CTGs.tpm       <- CTGs.tpm |> left_join(cancer.abbr) |> mutate(log2tpm = log2(median_tpm + 1), abbr = factor(abbr, levels = cancer.abbr$abbr))
# fig.1d         <- CTGs.tpm


CTGs.tpm      <- fig.1d
max.expr      <- max(CTGs.tpm$log2tpm)
proportion.df <- CTGs.tpm |> group_by(abbr) |> 
  summarise(proportion = mean(median_tpm > 1)) |> ungroup() |> 
  mutate(scaled_proportion = proportion * max.expr)


ggplot() +
  geom_boxplot(data = CTGs.tpm, aes(x = abbr, y = log2tpm), fill = "#F0F0F0", color = "#525252", outliers = FALSE) +
  geom_jitter(data = CTGs.tpm, aes(x = abbr, y = log2tpm), width = 0.1, alpha = 0.5, color = "#4292C6", size = 0.3) +
  geom_line(data = proportion.df, aes(x = abbr, y = scaled_proportion, group = 1), color = "#FEE08B", linewidth = 0.5, alpha = 0.8) +
  geom_point(data = proportion.df, aes(x = abbr, y = scaled_proportion), color = "#D53E4F", size = 1, shape = 17) +
  scale_y_continuous(name = expression(log[2]("TPM + 1")), breaks = c(0, 1, 2.5, 5, 7.5),
                     sec.axis = sec_axis(~ . / max.expr, name = "Proportion of CTGs with TPM >1", labels = scales::percent)) +
  labs(x = NULL, title = NULL) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#E7298A", linewidth = 0.3, alpha = 0.8)



## Fig. 1E-F ---------------------------------------------------------------------
library(readxl)
library(data.table)
library(GSVA)
library(BiocParallel)
library(foreach)


## CTG sets
hgnc         <- fread("./Datas/HGNC/non_alt_loci_set.txt")
hgnc.protein <- fread("./Datas/HGNC/gene_with_protein_product.txt")
CTdb         <- read_xlsx("./Proj_CTAs/input/CTGs/CTdatabase/CTAs.xlsx", skip = 2)
CTdb.ensembl <- hgnc |> filter(symbol %in% unique(CTdb$`Office symbol`), !grepl("not on reference assembly", location)) |> pull(ensembl_gene_id)


load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
CTGs.df      <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")


CTGs.Wang    <- read_xlsx("./Proj_CTAs/input/CTGs/Others/CTGs_Wang_NC_2016.xlsx", 2, skip = 1)
CTGs.Carter  <- read.csv("./Proj_CTAs/input/CTGs/Others/CTGs_Carter_JITC_2023.csv", header = TRUE)
CTGs.Carter  <- CTGs.Carter |> filter(!thymus_expression, tumor_expression)


CTGset       <- list(
  CTDatabase            = CTdb.ensembl,
  CTGs_Wang_NC_2016     = CTGs.Wang$ENname,
  CTGs_Carter_JITC_2023 = CTGs.Carter$ensembl_gene_id,
  candidates            = CTGs.df$ensembl,
  CTGs_us               = CTGs.df |> filter(symbol %in% pan.sample.CTGs) |> pull(ensembl)
)


## visualization
cancer.info <- TCGAbiolinks::getGDCprojects() |> filter(grepl("^TCGA", id)) |> select(id, cancer = name, abbr = tumor)
cancers     <- cancer.info$abbr
colors      <- setNames(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Set3"))(length(cancers)), cancers)


##### (1) E #####
# load("./Datas/UCSC_Xena/normal.tumor.RData")
# protein.genes <- gene.info |> filter(ensembl %in% hgnc.protein$ensembl_gene_id, id %in% rownames(log2.tpm)) |> pull(id)
# ssgsea        <- gsva(expr          = log2.tpm[protein.genes, ], 
#                       gset.idx.list = lapply(CTGset, \(x) gene.info |> filter(ensembl %in% x) |> pull(id)),
#                       method        = "ssgsea",
#                       abs.ranking   = FALSE,
#                       ssgsea.norm   = TRUE,
#                       BPPARAM       = MulticoreParam(40))
# 
# 
# auc.df <- foreach(c = names(sample.ls), .combine = 'rbind') %do% {
#   sample.info  <- sample.ls[[c]]
#   y            <- factor(sample.info$group, levels = c("normal", "tumor"))
#   ssgsea.score <- t(ssgsea)[sample.info$sample, , drop = FALSE]
#   auc          <- foreach(i = colnames(ssgsea.score), .combine = 'c') %do% {
#     pROC::roc(y, ssgsea.score[, i], direction = "<")$auc |> as.numeric()
#   }
#   auc.df       <- data.frame(abbr = c, auc = auc, CTGset = colnames(ssgsea.score))
# }
# auc.df |> group_by(CTGset) |> summarise(median(auc))  ## 0.91; 24 cancer types
# levels <- c("CTdatabase", "Wang et al. (2016)", "Carter et al. (2023)", "This Study", "Candidates")
# auc.df <- auc.df |>
#   mutate(CTGset = case_match(
#     CTGset,
#     "CTDatabase"            ~ "CTdatabase",
#     "CTGs_Wang_NC_2016"     ~ "Wang et al. (2016)",
#     "CTGs_Carter_JITC_2023" ~ "Carter et al. (2023)",
#     "candidates"            ~ "Candidates",
#     "CTGs_us"               ~ "This Study")) |>
#   mutate(CTGset = factor(CTGset, levels = levels))
# 
# 
# 
# ## TCGA/TAGET bulk RNA-seq data for cancer types in this study
# source("./Proj_CTAs/code/visualization.R")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")
# sample.df     <- foreach(i = names(sample.ls), .combine = 'rbind') %do% {sample.ls[[i]] |> mutate(cancer = i)}
# ## other TCGA cancers
# drop.cancers  <- c("GBM", "BRCA", "COAD", "HNSC", "KIRC", "LIHC", "LUAD", "OV", "PAAD", "PRAD", "SKCM")
# 
# 
# ssgsea.df <- data.frame(sample = colnames(ssgsea), ssgsea = ssgsea["CTGs_us", ], row.names = NULL)
# df        <- inner_join(sample.df, ssgsea.df, by = "sample")
# df$group  <- factor(df$group, levels = c("normal", "tumour"))
# 
# 
# auc.disc <- df |> 
#   group_by(cancer) |> 
#   group_modify(~ {
#     auc <- pROC::roc(.x$group, .x$ssgsea, direction = "<")$auc |> as.numeric()
#     data.frame(auc)
#   }) |> 
#   left_join(cancer.abbr) |>
#   mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) |>
#   ungroup() |>
#   select(abbr, auc) |>
#   mutate(cohort = "Discovery Cohort") |>
#   filter(abbr != "NB")
# auc.other <- auc.df |> filter(CTGset == "This Study") |> select(abbr, auc) |> filter(abbr %in% drop.cancers == FALSE) |> mutate(cohort = "TCGA Other Cancers")
# fig.1e    <- rbind(auc.disc, auc.other)


fig.1e |>
  ggplot(aes(x = cohort, y = auc)) +
  gghalves::geom_half_boxplot(aes(fill = cohort), side = "l", position = position_nudge(x = -0.1), outlier.alpha = 0,
                              fill = "gray80", alpha = 0.7, lwd = 0.2) +
  ggbeeswarm::geom_beeswarm(aes(fill = abbr), side = 1,
                            cex = 3.5, shape = 21, size = 2.5, stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = c(cancer.colours, colors[match(as.character(fig.1e$abbr[fig.1e$cohort == "TCGA Other Cancers"]), names(colors))])) +
  labs(x = NULL, y = "ROC-AUC", title = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 1, angle = 30),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", linewidth = 0.15) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))





##### (2) F #####
# tcga.purity     <- readRDS("./Datas/TCGA/pan-cancer/tcga.purity.rds")
# gdc.log2.fpkm   <- fread("./Datas/UCSC_Xena/GDC-PANCAN/GDC-PANCAN.htseq_fpkm.tsv") |> column_to_rownames("Ensembl_ID") |> as.matrix()  ## log2(x+1)
# gdc.sample.info <- fread("./Datas/UCSC_Xena/GDC-PANCAN/GDC-PANCAN.basic_phenotype.tsv")
# gene.info.v22   <- fread("./Datas/UCSC_Xena/GDC-PANCAN/gencode.v22.annotation.gene.probeMap") |> mutate(ensembl_gene_id = str_sub(id, end = 15))
# 
# 
# purity      <- tcga.purity |> 
#   select(sample = `Sample ID`, purity = CPE) |> 
#   mutate(purity = as.numeric(purity)) |> 
#   filter(!is.nan(purity))
# sample.info <- gdc.sample.info |> 
#   filter(program == "TCGA", sample %in% intersect(purity$sample, colnames(gdc.log2.fpkm))) |> 
#   mutate(abbr = str_extract(project_id, "(?<=TCGA-).*"))
# 
# 
# protein.genes <- gene.info.v22 |> filter(ensembl_gene_id %in% hgnc.protein$ensembl_gene_id, id %in% rownames(gdc.log2.fpkm)) |> pull(id)
# ssgsea        <- gsva(expr          = gdc.log2.fpkm[protein.genes, sample.info$sample],  
#                       gset.idx.list = lapply(CTGset, \(i) gene.info.v22 |> filter(ensembl_gene_id %in% i) |> pull(id)),
#                       method        = "ssgsea",
#                       abs.ranking   = FALSE,
#                       ssgsea.norm   = TRUE,
#                       BPPARAM       = MulticoreParam(40))
# 
# 
# df     <- t(ssgsea) |> as.data.frame() |> rownames_to_column("sample") |> inner_join(sample.info) |> inner_join(purity)
# rho.df <- foreach(i = rownames(ssgsea), .combine = 'rbind') %do% {
#   fit.ls <- split(df, df$abbr) |> lapply(\(x) {cor.test(x[[i]], x$purity, method = "spearman")})
#   rho    <- sapply(fit.ls, `[[`, "estimate")
#   p.val  <- sapply(fit.ls, `[[`, "p.value")
#   data.frame(rho, p.val, row.names = NULL) |> mutate(abbr = names(fit.ls), CTGset = i, .before = everything())
# }
# rho.df |> group_by(CTGset) |> summarise(median(rho))  ## 0.42; 21 cancer types
# levels <- c("CTdatabase", "Wang et al. (2016)", "Carter et al. (2023)", "This Study", "Candidates")
# rho.df <- rho.df |>
#   mutate(CTGset = case_match(
#     CTGset,
#     "CTDatabase"            ~ "CTdatabase",
#     "CTGs_Wang_NC_2016"     ~ "Wang et al. (2016)",
#     "CTGs_Carter_JITC_2023" ~ "Carter et al. (2023)",
#     "candidates"            ~ "Candidates",
#     "CTGs_us"               ~ "This Study")) |>
#   mutate(CTGset = factor(CTGset, levels = levels))
# 
# 
# ## TCGA/TAGET bulk RNA-seq data for cancer types in this study
# source("./Proj_CTAs/code/visualization.R")
# load("./Proj_CTAs/input/bulk RNA-seq datasets/tcga.target.gtex.RData")
# sample.df     <- foreach(i = names(sample.ls), .combine = 'rbind') %do% {sample.ls[[i]] |> mutate(cancer = i)}
# ## other TCGA cancers
# drop.cancers  <- c("GBM", "BRCA", "COAD", "HNSC", "KIRC", "LIHC", "LUAD", "OV", "PAAD", "PRAD", "SKCM")
# 
# 
# ssgsea.df        <- data.frame(sample = colnames(ssgsea), ssgsea = ssgsea["CTGs_us", ], row.names = NULL)
# ssgsea.purity.df <- inner_join(ssgsea.df, purity, by = "sample") |> mutate(sample_id = sample, sample = str_sub(sample_id, end = -2))
# dupli.samples    <- ssgsea.purity.df$sample[duplicated(ssgsea.purity.df$sample)]
# df               <- inner_join(sample.df, ssgsea.purity.df, by = "sample") |> filter(sample %in% dupli.samples == FALSE)
# df$group         <- factor(df$group, levels = c("normal", "tumour"))
# 
# 
# rho.disc <- df |> 
#   group_by(cancer) |> 
#   group_modify(~ {
#     fit <- cor.test(.x$ssgsea, .x$purity, method = "spearman")
#     data.frame(rho = fit$estimate, p.val = fit$p.value)
#   }) |> 
#   left_join(cancer.abbr) |>
#   mutate(abbr = factor(abbr, levels = cancer.abbr$abbr)) |>
#   ungroup() |>
#   select(abbr, rho, p.val) |>
#   mutate(cohort = "Discovery Cohort")
# rho.other <- rho.df |> filter(CTGset == "This Study") |> select(abbr, rho, p.val) |> filter(abbr %in% drop.cancers == FALSE) |> mutate(cohort = "TCGA Other Cancers")



load("../data/fig.1f.RData")  ## c("rho.disc", "rho.other", "df")
p1 <- rbind(rho.disc, rho.other) |>
  ggplot(aes(x = cohort, y = rho)) +
  gghalves::geom_half_boxplot(aes(fill = cohort), side = "l", position = position_nudge(x = -0.1), outlier.alpha = 0,
                              fill = "gray80", alpha = 0.7, lwd = 0.2) +
  ggbeeswarm::geom_beeswarm(aes(fill = abbr), side = 1,
                            cex = 3.5, shape = 21, size = 2.5, stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = c(cancer.colours, colors[match(rho.other$abbr, names(colors))])) +
  labs(x = NULL, y = expression(paste("Spearman's ", rho)), title = NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 1, angle = 30),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.25))
p2 <- df |>
  filter(cancer == "Breast Cancer (non_TNBC)") |>
  ggpubr::ggscatter(x = "ssgsea", y = "purity", color = "#2E5984", size = 1.5, alpha = 0.2,
                    add = "reg.line", add.params = list(color = "#C85200", fill = "gray"),
                    conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "spearman", size = 3.5)) + 
  labs(x = "CTG expression", y = "Tumor purity", title = "non_TNBC") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(linewidth = 0.1, color = "grey95"),
        axis.line = element_line(linewidth = 0.3, color = "black")) +
  coord_fixed(ratio = 1) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
p1 | p2

