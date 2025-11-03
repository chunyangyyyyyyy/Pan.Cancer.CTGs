library(tidyverse)
library(foreach)
library(doParallel)


source("./Proj_CTAs/code/visualization.R")
load("./Proj_CTAs/output/3. define CTGs/pan.sample.CTGs.RData")
CTGs.df <- readRDS("./Proj_CTAs/input/CTGs/CTGs.df.rds")



## Fig. 3A ---------------------------------------------------------------------
# load("./Proj_CTAs/output/4. characteristics of CTGs/genes.correlated.with.multi.CTGs.RData")
# rho.df <- foreach(i = cancer.abbr$cancer, .combine = 'rbind') %do% {
#   gene.CTG.cor.ls[[i]] %>%
#     mutate(gene_group = case_match(
#       gene,
#       pan.sample.CTGs.ls[[i]]          ~ "CTG",
#       setdiff(.$gene, pan.sample.CTGs) ~ "Non_CTG"
#     )) |>
#     filter(!is.na(gene_group), CTG %in% pan.sample.CTGs.ls[[i]], gene != CTG, !is.na(rho)) |>
#     group_by(CTG) |>
#     group_modify(~ {
#       CTG.median     <- .x |> filter(gene_group == "CTG")     |> pull(rho) |> median()
#       non.CTG.median <- .x |> filter(gene_group == "Non_CTG") |> pull(rho) |> median()
#       data.frame(cancer = i, rho = c(CTG.median, non.CTG.median), gene_group = c("CTG", "Non_CTG"))
#     }) |>
#     ungroup()
# }
# rho.df <- rho.df |> left_join(cancer.abbr) |> mutate(abbr = factor(abbr, levels = cancer.abbr$abbr))
# fig.3a <- rho.df


rho.df <- fig.3a
rho.df |>
  ggplot(aes(x = gene_group, y = rho)) +
  geom_boxplot(aes(color = gene_group), outliers = FALSE) +
  geom_point(aes(color = gene_group), size = 0.5, alpha = 0.4,
             position = position_jitterdodge(jitter.width = 0.8)) +
  ggpubr::stat_compare_means(group = "CTG", comparisons = list(c(c("CTG", "Non_CTG"))),
                             tip.length = 0, paired = TRUE, size = 2) +
  facet_wrap(~ abbr, ncol = 13) +
  scale_color_manual(values = c("Non_CTG" = "#4B8B3B", "CTG" = "#C0504D"), name = "Gene Group") +
  labs(x = "Discovery cohort", y = expression(paste("Spearman's ", rho))) +
  theme_minimal() +
  theme(strip.background    = element_rect(fill = "white"),
        strip.text          = element_text(face = "bold", size = 7),
        panel.grid.major.y  = element_line(color = "grey95", size = 0.2),
        legend.position     = "bottom",
        axis.text.x         = element_blank(),
        axis.ticks.x         = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))



## Fig. 3B ---------------------------------------------------------------------
library(igraph)
library(ggraph)


load("../data/fig.3b.RData")  ## c("cscore.ls", "gene.pair.ls")
colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(16)


## (1) the direction of co-expression for a gene-pair should be identical across all datasets
## (2) meidan(coex) > 0.25
i          <- "GBM"
df         <- gene.pair.ls[[i]] |> filter(coex > 0.25)
net        <- graph_from_data_frame(df, directed = FALSE)
cluster    <- cluster_louvain(net)
cluster.df <- data.frame(gene = cluster$names, cluster = as.character(cluster$membership), degree = degree(net))
hub.genes  <- cluster.df |> group_by(cluster) |> slice_max(order_by = degree, n = 2) |> pull(gene)
cluster.df <- cluster.df |> mutate(label = ifelse(gene %in% hub.genes, gene, NA), group = case_match(gene, unlist(Seurat::cc.genes.updated.2019) ~ "Cell Cycle (known)", .default = "Other"))


set.seed(1234)
p    <- ggraph(net, layout = 'igraph', algorithm = 'fr')
p1   <- p + 
  geom_edge_link(color = "lightgray", alpha = 0.5) + 
  geom_node_point(aes(fill = cluster.df$group, color = cluster.df$cluster), alpha = 0.7, size = 6, shape = 21, stroke = 0.4) + 
  scale_color_manual(values = colors, name = "Cluster") +
  scale_fill_manual(values = c("Cell Cycle (known)" = "firebrick", "Other" = "grey90"), name = NULL) +
  geom_node_text(label = cluster.df$label, repel = TRUE, color = "black", fontface = "italic", size = 6) +
  labs(title = i) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1)
p1



## Fig. 3C ---------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)


load("../data/fig.3c.RData")  ## "gene.CTG.cor.ls"
cancers                          <- names(pan.sample.CTGs.ls)
top2k.multiCTGs.cor.genes        <- foreach(i = cancers) %do% {
  rho.df <- gene.CTG.cor.ls[[i]] |> filter(CTG %in% pan.sample.CTGs.ls[[i]], gene %in% pan.sample.CTGs == FALSE) |> mutate(padj = p.adjust(pval, method = "BH"))
  genes  <- rho.df |> filter(padj < 0.01, rho > 0.5) |> count(gene) |> ungroup() |> slice_max(order_by = n, n = 2000) |> pull(gene)
}
names(top2k.multiCTGs.cor.genes) <- cancers


genesets   <- Filter(\(x) length(x) > 0, top2k.multiCTGs.cor.genes)
Entrez.ids <- lapply(genesets, \(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID)
x          <- compareCluster(Entrez.ids, fun = "enrichPathway", readable = TRUE)


p                      <- dotplot(x, font.size = 10, label_format = Inf, showCategory = 3) + 
  scale_size_continuous(range = c(3, 6)) + 
  labs(x = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
df                     <- p$data |> mutate(cancer = str_extract(Cluster, ".+(?=\\n)"))
levels(df$Description) <- str_wrap(levels(df$Description), width = 55)
cancer.omitted         <- setdiff(cancer.abbr$cancer, df$cancer)
df.add                 <- df[rep(1, length(cancer.omitted)), ] |> mutate(cancer = cancer.omitted, GeneRatio = NA, p.adjust = NA)
df                     <- rbind(df, df.add) |> full_join(cancer.abbr) |> mutate(Cluster = factor(abbr, levels = cancer.abbr$abbr))
p1                     <- p %+% df + labs(title = "Total CTGs") + theme(plot.title = element_text(hjust = 0.5))
p1



## Fig. 3D ---------------------------------------------------------------------
library(magrittr)
library(patchwork)


## (1) sc/snRNA-seq datasets
load("./Proj_CTAs/output/4. characteristics of CTGs/genes.correlated.with.multi.CTGs.RData")
x          <- x.all
p          <- clusterProfiler::dotplot(x, font.size = 10, label_format = Inf, showCategory = 3) 
df         <- p$data |> mutate(cancer = str_extract(Cluster, ".+(?=\\n)"))
## top2k genes involved in "Chromatin modifying enzymes" pathway
gene.ls    <- df |> filter(Description == "Chromatin modifying enzymes") %$% split(geneID, cancer) |> lapply(\(x) strsplit(x, "/")[[1]])
df1        <- data.frame(gene = unlist(gene.ls)) |> count(gene) |> mutate(rank = rank(-n, ties.method = "random"))
## (2) TCGA/TARGET bulk RNA-seq datasets
load("./Proj_CTAs/output/4. characteristics of CTGs/genes.correlated.with.multi.CTGs.bulk.RData")
x          <- x.all
p          <- clusterProfiler::dotplot(x, font.size = 10, label_format = Inf, showCategory = 3) 
df         <- p$data |> mutate(cancer = str_extract(Cluster, ".+(?=\\n)"))
## top2k genes involved in "Chromatin modifying enzymes" pathway
gene.ls    <- df |> filter(Description == "Chromatin modifying enzymes") %$% split(geneID, cancer) |> lapply(\(x) strsplit(x, "/")[[1]])
df2        <- data.frame(gene = unlist(gene.ls)) |> count(gene) |> mutate(rank = rank(-n, ties.method = "random"))


load("../data/fig.3d.RData")  ## c("df1", "df2")
epi.genes  <- readxl::read_xlsx("./Datas/Epifactors/v2.0/EpiGenes_main.xlsx")
set.seed(1234) 
df1.jitter <- df1 |> mutate(rank_j = rank + runif(n(), -0.3, 0.3), n_j = n + runif(n(), -0.15, 0.15))
df1.jitter <- df1.jitter |> left_join(epi.genes, by = c("gene" = "HGNC_symbol")) |> mutate(Target = ifelse(Target == "#" | is.na(Target), "unknown", Target))
set.seed(1234) 
df2.jitter <- df2 |> mutate(rank_j = rank + runif(n(), -0.3, 0.3), n_j = n + runif(n(), -0.15, 0.15))
df2.jitter <- df2.jitter |> left_join(epi.genes, by = c("gene" = "HGNC_symbol")) |> mutate(Target = ifelse(Target == "#" | is.na(Target), "unknown", Target))


top.regulators <- intersect(df1 |> slice_max(order_by = n, n = 10) |> pull(gene), 
                            df2 |> slice_max(order_by = n, n = 10) |> pull(gene))
p1 <- ggplot(df1.jitter, aes(x = rank_j, y = n_j)) +
  geom_point(aes(color = Target), alpha = 0.8, size = 1.5) +
  geom_point(data = filter(df1.jitter, gene %in% top.regulators),
             color = "firebrick", size = 3, shape = 21) +
  ggrepel::geom_text_repel(data = filter(df1.jitter, gene %in% top.regulators),
                           aes(label = gene), color = "firebrick", segment.color = "grey40", fontface = "italic",
                           box.padding = 0.6, min.segment.length = 0.2, max.overlaps = 100, 
                           force = 5, nudge_x = 0.15, nudge_y = 0.1) +
  scale_y_continuous(breaks = 1:10) +
  labs(x = "Gene Rank", y = "Number of Cancer Types", subtitle = "Discovery cohort") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p2 <- ggplot(df2.jitter, aes(x = rank_j, y = n_j)) +
  geom_point(aes(color = Target), alpha = 0.8, size = 1.5) +
  geom_point(data = filter(df2.jitter, gene %in% top.regulators),
             color = "firebrick", size = 3, shape = 21) +
  ggrepel::geom_text_repel(data = filter(df2.jitter, gene %in% top.regulators),
                           aes(label = gene), color = "firebrick", segment.color = "grey40", fontface = "italic",
                           box.padding = 0.6, min.segment.length = 0.2, max.overlaps = 100, 
                           force = 5, nudge_x = 0.15, nudge_y = 0.1) +
  scale_y_continuous(breaks = 1:10) +
  labs(x = "Gene Rank", y = "Number of Cancer Types", subtitle = "TCGA/TARGET bulk RNA-seq datasets") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))


p1 + p2 + plot_layout(ncol = 2, guides = "collect") + 
  plot_annotation(title = "Epigenetic Regulators Associated with CTG Activation", theme = theme(plot.title = element_text(hjust = 0.5)))



## Fig. 3E ---------------------------------------------------------------------
DE.res     <- read.table("./Proj_CTAs/input/bulk RNA-seq datasets/GSE90781_shYEATS2_H1299/GSE90781_H1299.shYEATS2-shCtrl.edgeR.txt.gz", header = TRUE)
geneList   <- setNames(DE.res$log2FC, DE.res$ID) |> sort(decreasing = TRUE)
term2gene1 <- CTGs.df |> filter(symbol %in% pan.sample.CTGs) |> mutate(term = "Total CTGs (n = 407)") |> dplyr::select(term, gene = ensembl)
term2gene2 <- CTGs.df |> filter(symbol %in% pan.sample.CTGs.ls$`Lung Cancer (LUAD)`) |> mutate(term = "LUAD CTGs (n = 190)") |> dplyr::select(term, gene = ensembl)
term2gene  <- rbind(term2gene1, term2gene2)


## Warning: The order of those tied genes will be arbitrary, which may produce unexpected results.
set.seed(123)
gsea <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene, pvalueCutoff = 1)
enrichplot::gseaplot2(gsea, geneSetID = 1:2, pvalue_table = TRUE,
                      color = c("#E495A5", "#86B875"), ES_geom = "line", title = "H1299: YEATS2-KD vs NT-Ctrl")



## Fig. 3F ---------------------------------------------------------------------
# cl <- makeCluster(20)
# registerDoParallel(cl)
# 
# 
# sample.info         <- readRDS("./Proj_CTAs/input/scATAC-seq datasets/sample.info.rds")
# log.tp10k.ls        <- foreach(i = sample.info$rds_file, .packages = c('Seurat', 'Signac')) %dopar% {
#   tryCatch(
#     expr  = {
#       rds.path <- file.path("./Proj_CTAs/input/scATAC-seq datasets/processed", i)
#       seu      <- readRDS(rds.path)
#       stopifnot(DefaultAssay(seu) == "RNA")
#       seu.mal  <- subset(seu, cell_type == "Tumor")
#       seu.mal  <- NormalizeData(seu.mal)
#       seu.mal[["RNA"]]@data
#     },
#     error = function(e) e$message
#   )
# }
# names(log.tp10k.ls) <- sample.info$Piece_ID
# idx                 <- sapply(log.tp10k.ls, \(x) class(x) == "dgCMatrix")
# log.tp10k.ls        <- log.tp10k.ls[idx]
# sapply(log.tp10k.ls, \(x) identical(rownames(x), rownames(log.tp10k.ls[[1]]))) |> all() |> stopifnot()
# 
# 
# cancers               <- unique(sample.info$cancer)
# cancer.expr.ls        <- foreach(c = cancers) %do% {
#   s    <- sample.info |> filter(cancer == c) |> pull(Piece_ID) |> intersect(names(log.tp10k.ls))
#   ls   <- log.tp10k.ls[s]
#   ls   <- foreach(i = s) %do% {
#     x           <- ls[[i]]
#     colnames(x) <- paste(i, colnames(x), sep = "_")
#     x
#   }
#   expr <- Reduce(cbind, ls)
# }
# names(cancer.expr.ls) <- cancers
# 
# 
# regulators    <- c("ARID4B", "KDM3B", "PBRM1", "SETD2", "SUZ12", "YEATS2")
# rho.ls        <- foreach(c = cancers) %do% {
#   expr          <- cancer.expr.ls[[c]] |> as.matrix() |> t()
#   rho.ls        <- foreach(i = regulators) %do% {
#     x    <- expr[, i]
#     rho  <- cor(expr, x, method = "spearman")
#     pval <- correlation::cor_to_p(rho, nrow(expr), method = "spearman")$p
#     
#     rho.df <- data.frame(gene = rownames(rho), rho = rho[, 1], p.val = pval[, 1], row.names = NULL)
#     rho.df <- rho.df |> 
#       filter(gene != i) |> 
#       mutate(p.adj = p.adjust(p.val, method = "fdr"),
#              rho   = ifelse(p.adj >= 0.05, 0, rho))
#   }
#   names(rho.ls) <- regulators
#   
#   rho.ls
# }
# names(rho.ls) <- cancers
# 
# 
# cancers   <- names(cancer.expr.ls)
# term2gene <- data.frame(term = "Total CTGs (n = 407)", gene = pan.sample.CTGs)
# gsea.df   <- foreach(c = cancers, .combine = 'rbind') %do% {
#   ls      <- rho.ls[[c]]
#   gsea.df <- foreach(i = names(ls), .combine = 'rbind') %do% {
#     rho.df   <- ls[[i]]
#     geneList <- setNames(rho.df$rho, rho.df$gene) |> sort(decreasing = TRUE)
#     set.seed(1234)
#     gsea     <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene, pvalueCutoff = 1)
#     
#     data.frame(regulator = i, NES = gsea@result$NES, p.val = gsea@result$pvalue)
#   }
#   gsea.df |> mutate(cancer = c)
# }
# fig.3f    <- gsea.df


gsea.df <- fig.3f
p.sig   <- 0.05
df      <- gsea.df |> 
  filter(cancer != "MM") |>
  mutate(sig_star   = case_when(
    p.val <= 0.0001 ~ "****",
    p.val <= 0.001 ~ "***",
    p.val <= 0.01 ~ "**",
    p.val <= 0.05 ~ "*"))
p       <- ggplot(df, aes(x = cancer, y = regulator, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.2, width = 0.9, height = 0.9) +
  geom_text(aes(label = sig_star),
            color = "black", size = 2.5, vjust = 0.7, fontface = "bold") +
  scale_fill_gradientn(colours = c("#F7F7F7", "#D9D9D9", "#F7CACA", "#F1898F", "#D94A6E", "#A62F5E"),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 0), 
        axis.text.y = element_text(face = "italic"), 
        panel.grid = element_blank())
p