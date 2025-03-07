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

# load packages
invisible(lapply(packages, library, character.only = T))

source("scripts/custom_functions.R")




# Plotting Options
group_colors <- c(
    "0d_SPF" = "#eafaf1",
    "0d_GF" = "#fef9e7",
    "3d_SPF_sal" = "#aed6f1",
    "3d_SPF_blm" = "#2e86c1",
    "3d_GF_sal" = "#f5b7b1",
    "3d_GF_blm" = "#e74c3c",
    "21d_SPF_sal" = "#a9dfbf",
    "21d_SPF_blm" = "#229954",
    "21d_GF_sal" = "#f0b27a",
    "21d_GF_blm" = "#af601a"
)

group_colors.transparent <- sapply(group_colors,
                                   function(x) grDevices::adjustcolor(x, alpha.f = .4))

group_shape <- c("No Treatment" = 22,  "Saline" = 21, "Bleomycin" = 24)


make_annot <- function(df, which, show_legend) {
    ha <- ComplexHeatmap::HeatmapAnnotation(df = df[, c("Duration", "MouseType", "Treatment")],
                                      col = list(Duration = c("0 Day" = "#eafaf1", "3 Days" = "#d2b4de", "21 Days" = "#f4d03f"),
                                                 MouseType = c("SPF" = "#58d68d", "GF" = "#dc7633"),
                                                 Treatment = c("No Treatment" = "#d5dbdb", "Saline" = "#5dade2", "Bleomycin"  = "#ec7063")),
                                      which = which,
                                      border = TRUE,
                                      show_annotation_name = FALSE,
                                      show_legend = show_legend)
    
    return(ha)
}


#@@@
# 0. save biomart query
#@@@
CacheBiomaRt()



#@@@
# 1. create DESeq object using salmon output
#@@@

# set base file path 
salmon_output_root <- "data/salmon_output"

# list up sample names
sample_names.orig <- list.files(salmon_output_root)

# change name order (no need to make it as factor)
sample_names.orig <- c(sample_names.orig[5:8],    #  0 day,  SPF, no treatment
                       sample_names.orig[1:4],    #  0 day,  GF,  no treatment
                       sample_names.orig[37:40],  #  3 days, SPF, saline
                       sample_names.orig[33:36],  #  3 days, SPF, bleomycin
                       sample_names.orig[29:32],  #  3 days, GF,  saline
                       sample_names.orig[25:28],  #  3 days, GF,  bleomycin
                       sample_names.orig[21:24],  # 21 days, SPF, saline
                       sample_names.orig[17:20],  # 21 days, SPF, bleomycin
                       sample_names.orig[13:16],  # 21 days, GF,  saline
                       sample_names.orig[9:12])   # 21 days, GF,  bleomycin

salmon_output_path <- file.path(salmon_output_root, sample_names.orig, "quant.sf")
names(salmon_output_path) <- sample_names.orig

# prepare tx2gene
ah <- AnnotationHub::AnnotationHub()
ahDb <- AnnotationHub::query(x = ah, pattern = c("ensDb", "Mus musculus"))

ah_latest_rec <- ahDb %>%
    mcols() %>%  # retrieve dataframe containing metadata columns
    rownames() %>%
    tail(n = 1)

ensDb <- ahDb[[ah_latest_rec]]

txData <- ensembldb::transcripts(x = ensDb, return.type = "DataFrame")
tx2gene <- txData[, c("tx_id", "gene_id")]

# create txi object
txi <- tximport::tximport(files = salmon_output_path, 
                          type = "salmon", 
                          tx2gene = tx2gene, 
                          ignoreTxVersion = T)

# prepare metadata
coldata <- data.frame(Sample = sample_names.orig)
coldata <- coldata %>%
    mutate(Group = sapply(sample_names.orig, function(x) {substr(x, 1, nchar(x)-1)}),
           Duration = case_when(grepl("0d", Sample) ~ "0 Day",
                                grepl("3d", Sample) ~ "3 Days",
                                grepl("21d", Sample) ~ "21 Days"),
           MouseType = case_when(grepl("SPF", Sample) ~ "SPF",
                                 grepl("GF", Sample) ~ "GF"),
           Treatment = case_when(grepl("sal", Sample) ~ "Saline",
                                 grepl("blm", Sample) ~ "Bleomycin",
                                 TRUE ~ "No Treatment")) %>%
    tibble::column_to_rownames(var = "Sample")

coldata$Group <- factor(coldata$Group, 
                        levels = c("0d_SPF", "0d_GF", 
                                   "3d_SPF_sal", "3d_SPF_blm", 
                                   "3d_GF_sal", "3d_GF_blm", 
                                   "21d_SPF_sal","21d_SPF_blm", 
                                   "21d_GF_sal", "21d_GF_blm"))

coldata$Duration <- factor(coldata$Duration,
                            levels = c("0 Day", "3 Days", "21 Days"))

coldata$MouseType <- factor(coldata$MouseType,
                            levels = c("SPF", "GF"))

coldata$Treatment <- factor(coldata$Treatment,
                            levels = c("No Treatment", "Saline", "Bleomycin"))

#@@@
# 2. basic sample QC & EDA
#@@@

# create DESeq object
dds <- DESeq2::DESeqDataSetFromTximport(txi = txi, 
                                        colData = coldata, 
                                        design = ~Group)

dds.clean <- dds <- dds[rowSums(BiocGenerics::counts(dds) >= 10) > 4, ]

ddsX <- DESeq2::DESeq(dds)

vst <- DESeq2::vst(ddsX, blind = F)


RLEplot(vst, group_colors = group_colors, box.only = TRUE)

PCAplot(vst,
        ntop = 500,
        grp_col = group_colors,
        pt_sz = 8, lab_sz = 4, glob_txt_sz = 16, leg_pos = "none") +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 2)

CorHeatmap(vst, 
           top_annot = make_annot(colData(dds), which = "column", show_legend = TRUE),
           right_annot = make_annot(colData(dds), which = "row", show_legend = FALSE), 
           nr_col_split = 7, nr_row_split = 7)

#------------------------------------------------------------------------------#
dds.SPF <- dds.clean[, grepl("SPF", colnames(dds.clean)) & !(colnames(dds.clean) %in% c("21d_SPF_sal1", "0d_SPF1", "21d_SPF_blm1", "3d_SPF_blm5"))]
dds.SPF$Group <- droplevels(dds.SPF$Group)
ddsX.SPF <- DESeq2::DESeq(dds.SPF)
vst.SPF <- DESeq2::vst(ddsX.SPF)

PCAplot(vst.SPF,
        ntop = 500,
        grp_col = group_colors,
        pt_sz = 8, lab_sz = 4, glob_txt_sz = 16, leg_pos = "right") +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 2)

CorHeatmap(vst.SPF, 
           top_annot = make_annot(colData(dds.SPF), which = "column", show_legend = TRUE),
           right_annot = make_annot(colData(dds.SPF), which = "row", show_legend = FALSE), 
           nr_col_split = 5, nr_row_split = 5)

#------------------------------------------------------------------------------#
dds.GF <- dds.clean[, grepl("GF", colnames(dds.clean)) & !(colnames(dds.clean) %in% c("21d_GF_sal5", "21d_GF_blm2", "3d_GF_sal1", "0d_GF1", "3d_GF_blm4"))]
dds.GF <- dds[, grepl("GF", colnames(dds)) & !(colnames(dds) %in% c("21d_GF_sal5"))]
dds.GF$Group <- droplevels(dds.GF$Group)
ddsX.GF <- DESeq2::DESeq(dds.GF)
vst.GF <- DESeq2::vst(ddsX.GF)

PCAplot(vst.GF,
        ntop = 500,
        grp_col = group_colors,
        pt_sz = 8, lab_sz = 4, glob_txt_sz = 16, leg_pos = "right") +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 2)

CorHeatmap(vst.GF, 
           top_annot = make_annot(colData(dds.GF), which = "column", show_legend = TRUE),
           right_annot = make_annot(colData(dds.GF), which = "row", show_legend = FALSE),
           nr_col_split = 5, nr_row_split = 5)

#------------------------------------------------------------------------------#
outliers <- c("0d_SPF1", "0d_GF1", "3d_SPF_blm5", "3d_GF_sal1", "3d_GF_blm4", "21d_SPF_sal1", "21d_SPF_blm1", "21d_GF_sal5", "21d_GF_blm2")

dds.clean <- dds.clean[, !(colnames(dds.clean) %in% outliers)]
ddsX.clean <- DESeq2::DESeq(dds.clean)
vst.clean <- DESeq2::vst(ddsX.clean, blind = F)

PCAplot(vst.clean,
        ntop = 500,
        grp_col = group_colors,
        pt_sz = 8, lab_sz = 4, glob_txt_sz = 16, leg_pos = "right") +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 2)


CorHeatmap(vst.clean, 
           top_annot = make_annot(colData(dds.clean), which = "column", show_legend = TRUE),
           right_annot = make_annot(colData(dds.clean), which = "row", show_legend = FALSE), 
           nr_col_split = 8, nr_row_split = 8)

# save outlier removed & ID translated read count data
cnt <- Ens2Sym(counts(dds.clean) %>% as.data.frame() %>% tibble::rownames_to_column(var = "ENSG"), "ENSG") %>%
    dplyr::select(Gene, ENSG, dplyr::everything())

writexl::write_xlsx(cnt, "results/raw_read_counts.xlsx")

# save design matrix of outlier removed data
writexl::write_xlsx(colData(dds.clean) %>% as.data.frame() %>% tibble::rownames_to_column(var = "Sample"), "results/design_matrix2.xlsx")

#------------------------------------------------------------------------------#




