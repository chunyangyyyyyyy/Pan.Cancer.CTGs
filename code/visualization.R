cancer.abbr <- tibble::tribble(
  ~ cancer,                                ~ abbr,
  "Brain Cancer (glioblastoma)",           "GBM",      
  "Breast Cancer (TNBC)",                  "TNBC",     
  "Breast Cancer (non_TNBC)",              "non_TNBC", 
  "Colorectal Cancer",                     "CRC",      
  "Head and Neck Cancer (HNSCC)",          "HNSCC",    
  "Kidney Cancer (clear cell carcinoma)",  "ccRCC",  
  "Liver Cancer (HCC)",                    "HCC", 
  "Lung Cancer (LUAD)",                    "LUAD",  
  "Neuroendocrine Cancer (neuroblastoma)", "NB",
  "Ovarian Cancer (HGSOC)",                "HGSOC",
  "Pancreatic Cancer (PDAC)",              "PDAC",  
  "Prostate Cancer",                       "PRAD",  
  "Skin Cancer (melanoma)",                "SKCM",
)


cancer.colours        <- c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFED6F", "#999999")
names(cancer.colours) <- cancer.abbr$abbr
