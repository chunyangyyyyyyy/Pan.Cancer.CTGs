library(Seurat)


calculate.n.genes.control <- function(seu, feature) {
  CreateRandomGenes <- function (seu, feature, nbin = 24, seed = 1) {
    ## the code was inspired from `AddModuleScore {Seurat}`
    
    set.seed(seed)
    
    data.avg <- Matrix::rowMeans(x = seu[["RNA"]]@data)
    data.avg <- data.avg[order(data.avg)]
    data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                                    n = nbin, labels = FALSE, right = FALSE)
    names(x = data.cut) <- names(x = data.avg)
    
    feature   <- intersect(feature, names(data.cut))
    ref.genes <- setdiff(names(data.cut), feature)
    
    nfeature.each.bin  <- as.list(table(data.cut[feature]))
    ref.genes.each.bin <- split(names(data.cut[ref.genes]), data.cut[ref.genes])
    
    random.genes <- lapply(names(nfeature.each.bin), \(i) {
      
      sample(x       = ref.genes.each.bin[[i]], 
             size    = nfeature.each.bin[[i]], 
             replace = FALSE)
    }
    )
    random.genes <- as.character(unlist(random.genes))
    
    return(random.genes)
  }
  
  
  set.seed(1234)
  seeds           <- sample(1000, 100)
  random.genesets <- lapply(seeds, \(i) CreateRandomGenes(seu, feature, seed = i))
  
  
  n.genes.control <- sapply(random.genesets, \(genes) colSums(seu[["RNA"]]@data[genes, , drop = FALSE] > 0))
  n.genes.control <- apply(n.genes.control, 1, median)
  
  
  return(n.genes.control)
}
