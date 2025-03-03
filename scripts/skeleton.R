# set working directory
if (strsplit(getwd(), '/')[[1]][length(strsplit(getwd(), '/')[[1]])] == 'scripts') {
    setwd("../")
}

# package handling
packages <- c("tidyverse", "reshape2", "readxl", "writexl",
              "AnnotationHub", "ensembldb", "biomaRt", "fst", "msigdbr",
              "sva", "tximport", "DESeq2", "apeglm", "ashr",
              "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db",
              "GSVA", "decoupleR", "OmnipathR", "clusterProfiler",
              "EnhancedVolcano", "circlize", "ggrepel", "ComplexHeatmap")

# install pacakges
BiocManager::install(packages)

# re-init renv
renv::init()  # select re-initialization