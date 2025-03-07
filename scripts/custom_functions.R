#===============================================#
#   CacheBiomaRt : cache biomart query result   #
#===============================================#
CacheBiomaRt <- function() {
    # ensembl_gene_id, external_gene_name, entrezgene_id, ensembl_transcript_id : "feature_page" page
    # ensembl_gene_id, hsapiens_homolog_ensembl_gene, hsapiens_homolog_associated_gene_name : "homologs" page
    # notice that "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id" are also in a homologs page
    queries.feature_page <- c("ensembl_gene_id",
                              "external_gene_name",
                              "entrezgene_id",
                              "ensembl_transcript_id")
    
    queries.homologs <- c("ensembl_gene_id", 
                          "hsapiens_homolog_ensembl_gene",
                          "hsapiens_homolog_associated_gene_name")
    
    bm.mart <- biomaRt::useEnsembl(biomart = "genes",
                                   dataset = "mmusculus_gene_ensembl")
    
    # different attribute pages can't be queried with a single getBM()
    bm.query.feature_page <- biomaRt::getBM(attributes = queries.feature_page,
                                            mart = bm.mart)
    
    bm.query.homologs <- biomaRt::getBM(attributes = queries.homologs,
                                        mart = bm.mart)
    
    # do not handle NA or empty values here. retain full data 
    bm.query <- merge(bm.query.feature_page, bm.query.homologs, by = "ensembl_gene_id", all = T)
    
    # always remember that bm.query will contain empty values, NA values, duplicated names
    assign("biomart.cache", bm.query, envir = .GlobalEnv)
    fst::write.fst(bm.query, "data/biomart.cache")
}



#=============================================#
#   Ens2Sym : Map Ensembl ID to Gene Symbol   #
#=============================================#
#======================================================================================================#
#   https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html   #
#======================================================================================================#
Ens2Sym <- function(df, ENSG_col) {
    
    if (!exists("biomart.cache")) {
        if (!file.exists("data/biomart.cache")) {
            CacheBiomaRt()
        } else {
            biomart.cache <- fst::read.fst("data/biomart.cache")
        }
    }
    
    BM <- biomart.cache %>%
        dplyr::filter(ensembl_gene_id %in% df[[ENSG_col]]) %>%
        dplyr::select(ensembl_gene_id, external_gene_name) %>%
        dplyr::distinct()
    
    colnames(BM) <- c("EnsemblID", "Gene")
    
    # remove redundant genes
    # unmapped genes --> ''
    # remove pseudo-genes, unannotated genes --> '-ps', 'Gm', 'Rik'
    # remove anti-sense genes --> '-as'
    # remove opposite-strand genes --> '[0-9]os'
    # remove ribosomal genes --> 'Rps', 'Rpl'
    # remove mitochondrial genes --> 'mt-'
    # remove ORF clone --> 'AB1234'
    df <- df %>% 
        dplyr::left_join(BM, by = c("ENSG" = "EnsemblID")) %>%
        dplyr::filter(
            (Gene != '') & !grepl("[0-9]Rik", Gene) & !grepl("Gm[0-9]", Gene) & !grepl("-ps", Gene) & !grepl("-as", Gene) & !grepl("Rp[ls]", Gene) & !grepl("mt-", Gene) & !grepl("[0-9]os", Gene) & !grepl("[A-Z][A-Z][0-9]", Gene)
        ) %>% drop_na()
    
    return(df)
}



#=============================#
#   RLEplot : Draw RLE plot   #
#=============================#
RLEplot <- function(X_vst, group_colors, raw.mat = FALSE, font_size = 15, box.only = F) {
    if (raw.mat) {
        # in this case, X_vst is raw matrix count even if the name is X_vst
        dat <- log2(X_vst + 1)
    } else {
        dat <- log2(SummarizedExperiment::assay(X_vst) + 1)
    }
    
    medians_per_gene <- Biobase::rowMedians(dat)  # median values per gene
    
    # make median of each gene equals to 0
    RLE <- dat - medians_per_gene  # (gene) x (sample)
    
    # make wide-format to long-format
    RLE_long <- reshape2::melt(RLE) %>%
        dplyr::mutate(Group = dplyr::case_when(grepl("0d_SPF", Var2) ~ "0d_SPF",
                                               grepl("0d_GF", Var2) ~ "0d_GF",
                                               grepl("3d_SPF_sal", Var2) ~ "3d_SPF_sal",
                                               grepl("3d_SPF_blm", Var2) ~ "3d_SPF_blm",
                                               grepl("3d_GF_sal", Var2) ~ "3d_GF_sal",
                                               grepl("3d_GF_blm", Var2) ~ "3d_GF_blm",
                                               grepl("21d_SPF_sal", Var2) ~ "21d_SPF_sal",
                                               grepl("21d_SPF_blm", Var2) ~ "21d_SPF_blm",
                                               grepl("21d_GF_sal", Var2) ~ "21d_GF_sal",
                                               grepl("21d_GF_blm", Var2) ~ "21d_GF_blm"))
    
    RLE_long$Group <- factor(RLE_long$Group,
                             levels = c("0d_SPF", "0d_GF", 
                                        "3d_SPF_sal", "3d_SPF_blm", 
                                        "3d_GF_sal", "3d_GF_blm", 
                                        "21d_SPF_sal", "21d_SPF_blm", 
                                        "21d_GF_sal", "21d_GF_blm"))
    
    colnames(RLE_long) <- c("ENSG", "Sample", "RLE", "Group")
    test_a <<- RLE_long
    if (box.only) {
        RLE_plot <- ggplot(data = RLE_long, aes(x = Sample, y = RLE)) + 
            geom_hline(yintercept = 0, color = "red", linewidth = 1) +
            geom_boxplot(aes(fill = Group), outlier.shape = NA, coef = 0) +
            scale_fill_manual(values = group_colors) +
            theme_bw() +
            theme(text = element_text(size = font_size),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            xlab("Samples") + ylab("Relative Log Expression")
        
        # calculate sample-wise box limits
        upper_bound <- max(apply(RLE, MARGIN = 2, function(x) quantile(x, 0.75)))
        lower_bound <- min(apply(RLE, MARGIN = 2, function(x) quantile(x, 0.25)))
        
        RLE_plot <- RLE_plot + 
            coord_cartesian(ylim = c(lower_bound - 0.001, upper_bound + 0.001))
        
    }else {
        RLE_plot <- ggplot(data = RLE_long, aes(x = Sample, y = RLE)) + 
            geom_hline(yintercept = 0, color = "red", linewidth = 1) +
            geom_boxplot(aes(fill = Group)) +
            scale_fill_manual(values = group_colors) +
            theme_bw() + 
            theme(text = element_text(size = font_size), 
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            xlab("Samples") + ylab("Relative Log Expression")
    }
    
    return(RLE_plot)
}



#===============================#
#   PCAplot : Plot PCA result   #
#===============================#
PCAplot <- function(X_vst, ntop = 500, genes = NULL, metadata = NULL,
                    grp = "Group", grp_col, ext_grp = NULL, grp_shape = NULL,
                    pt_sz = 5, lab_sz = 3, glob_txt_sz = 10, leg_pos = "right") {
    if (is(X_vst, "DESeqTransform")){
        dat <- SummarizedExperiment::assay(X_vst)
        group <- X_vst[[grp]]
        names <- colnames(X_vst)
        metadata <- colData(X_vst)
    } else {
        if (is.null(metadata))  stop("You must give metadata if you want to use custom dataframe!")
        group <- metadata[[grp]]
        names <- rownames(metadata)
        
        if (is.matrix(X_vst)) {
            dat <- X_vst
        } else if (is.data.frame(X_vst)) {
            dat <- as.matrix(X_vst)
        }
        
        dat <- log2(dat+1)
    }
    
    if (is.null(genes)) {
        rv <- rowVars(dat)
        hvg <- order(rv, decreasing = T)[seq_len(min(ntop, length(rv)))]
        dat <- t(dat[hvg, ])  # [matrix] (sample) x (features)
    } else {
        genes <-biomart.cache %>% 
            dplyr::filter(external_gene_name %in% genes) %>% 
            dplyr::pull(ensembl_gene_id) %>%
            unique()
        
        dat <- t(dat[rownames(dat) %in% genes, ])  # [matrix] (sample) x (features)
    }
    
    if (!all(apply(dat, MARGIN = 2, var) != 0)) {
        dat <- dat[, which(apply(dat, MARGIN = 2, var) != 0)]
    }
    
    pca <- stats::prcomp(dat)
    
    percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
    
    if (is.null(ext_grp)) {
        pca.df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                             Group = group)
        
        rownames(pca.df) <- names
        
        plot <- ggplot(pca.df, aes(x = PC1, y = PC2)) +
            geom_point(size = pt_sz, aes(fill = Group), shape = 21) + 
            scale_fill_manual(values = grp_col) + 
            xlab(paste0("PC1 : ", percentVar[1], '%')) +
            ylab(paste0("PC2 : ", percentVar[2], '%')) +
            geom_text_repel(size = lab_sz, label = rownames(pca.df)) +
            theme_minimal() + 
            theme(text = element_text(size = glob_txt_sz), legend.position = leg_pos)
    } else {
        if (is.null(grp_shape)) {
            stop("you should give group shape!")
        }
        
        pca.df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                             Group = group,
                             ExtGrp = metadata[ext_grp][, 1])
        
        rownames(pca.df) <- names
        
        plot <- ggplot(pca.df, aes(x = PC1, y = PC2)) +
            
            # you should use shape 21 ~ 25 to use both shape and fill
            # https://stackoverflow.com/questions/15965870/fill-and-border-colour-in-geom-point-scale-colour-manual-in-ggplot
            geom_point(size = pt_sz, aes(fill = Group, shape = ExtGrp)) + 
            
            scale_shape_manual(values = grp_shape) +
            scale_fill_manual(values = grp_col) + 
            xlab(paste0("PC1 : ", percentVar[1], '%')) +
            ylab(paste0("PC2 : ", percentVar[2], '%')) +
            geom_text_repel(size = lab_sz, label = rownames(pca.df)) +
            theme_minimal() + 
            theme(text = element_text(size = glob_txt_sz), legend.position = leg_pos)
    }
    
    return(plot)
}



#=======================================================#
#   CorHeatmap : Plot sample-wise correlation heatmap   #
#=======================================================#
CorHeatmap <- function(X_vst, 
                       nr_col_split = NULL,
                       nr_row_split = NULL,
                       top_annot = NULL, bottom_annot = NULL,
                       left_annot = NULL, right_annot = NULL) {
    
    cor_mat <- stats::cor(SummarizedExperiment::assay(X_vst))
    
    heatmap <- ComplexHeatmap::Heatmap(matrix = cor_mat,
                                       rect_gp = gpar(col = "white", lwd = 0.5),
                                       row_split = nr_row_split,
                                       column_split = nr_col_split,
                                       row_title = NULL,
                                       column_title = NULL,
                                       top_annotation = top_annot,
                                       bottom_annotation = bottom_annot,
                                       left_annotation = left_annot,
                                       right_annotation = right_annot,
                                       show_heatmap_legend = T,
                                       heatmap_legend_param = list(
                                           title = "Corr"
                                       ))
    
    return(heatmap)
}



#================================================================#
#   MarkersExprsHeatmap : Plot marker genes expression heatmap   #
#================================================================#
MarkersExprsHeatmap <- function(X_vst,
                                ENSG_Markers,
                                scale = T,
                                legend_name = "Z-score",
                                custom_col = NULL,
                                show_row_names = F,
                                cluster_col = T,
                                cluster_row = T,
                                col_annot = NULL,
                                row_annot = NULL,
                                nr_row_split = NULL,
                                nr_col_split = NULL) {
    
    dat <- SummarizedExperiment::assay(X_vst)  # [matrix] (features) x (samples)
    
    ENSG_Markers <- ENSG_Markers[ENSG_Markers %in% rownames(dat)]
    markers <- biomart.cache %>% 
        dplyr::filter(ensembl_gene_id %in% ENSG_Markers) %>%
        dplyr::distinct(ensembl_gene_id, .keep_all = T) %>%
        dplyr::pull(external_gene_name)
    
    dat <- dat[ENSG_Markers, ]
    rownames(dat) <- markers
    
    test.dat <<- dat
    
    if (!scale) {
        norm_dat <- dat
    } else {
        norm_dat <- dat %>%
            t() %>%  # (samples) x (features)
            scale() %>%  # center and scale each column (gene) (row scaling of original data)
            t()  # (features) x (samples)
    }
    
    heatmap <- ComplexHeatmap::Heatmap(matrix = norm_dat,
                                       rect_gp = gpar(col = "white", lwd = 0.5),
                                       col = custom_col,
                                       cluster_columns = cluster_col,
                                       cluster_rows = cluster_row,
                                       show_row_names = show_row_names,
                                       column_names_rot = 45,  # rotating sample names 45 degrees
                                       row_split = nr_row_split,
                                       row_title = NULL,
                                       column_split = nr_col_split,
                                       column_title = NULL,
                                       border = T,
                                       bottom_annotation = col_annot,
                                       right_annotation = row_annot,
                                       use_raster = F,
                                       name = legend_name)
    
    return(heatmap)
}



#=================================================================#
#   runCollecTRI : run Multi-Linear Model on CollecTRI database   #
#=================================================================#
runCollecTRI <- function(X_vst) {
    ctri.db <- decoupleR::get_collectri(organism = "mouse")
    
    if (!exists("biomart.cache")) {
        if (!file.exists("data/biomart.cache")) {
            CacheBiomaRt()
        } else {
            biomart.cache <- fst::read.fst("data/biomart.cache")
        }
    }
    
    ctri.ens2sym <- biomart.cache %>%
        dplyr::filter(ensembl_gene_id %in% rownames(X_vst)) %>%
        dplyr::select(ensembl_gene_id, external_gene_name) %>%
        dplyr::distinct()
    
    ctri.dat <- SummarizedExperiment::assay(X_vst) %>%
        as.data.frame() %>%
        dplyr::mutate_if(~ any(is.na(.x)), ~ if_else(is.na(.x), 0, .x)) %>%
        tibble::rownames_to_column(var = "EnsemblID") %>%
        dplyr::inner_join(ctri.ens2sym, by = c("EnsemblID" = "ensembl_gene_id")) %>%
        dplyr::filter(external_gene_name != '') %>%
        dplyr::mutate(gene_means = rowMeans(dplyr::across(-c("EnsemblID", "external_gene_name")))) %>%
        dplyr::arrange(dplyr::desc(gene_means)) %>%
        dplyr::distinct(external_gene_name, .keep_all = T) %>%
        dplyr::select(-EnsemblID, -gene_means) %>%
        tibble::column_to_rownames(var = "external_gene_name") %>%
        as.matrix()
    
    ctri.result <- decoupleR::run_ulm(mat = ctri.dat,
                                      net = ctri.db,
                                      .source = "source",
                                      .target = "target",
                                      .mor = "mor",
                                      minsize = 5)
    
    #ctri.result.visual_form <- ctri.result %>%
    
    return(ctri.result)
}



#============================================#
#   runGSVA : run GSVA on hallmark geneset   #
#============================================#
#=======================================================================================================#
#   https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html                    #
#   https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_03_gsva.html   #
#=======================================================================================================#
runGSVA <- function(X_vst, target_gs = NULL) {
    if (is.null(target_gs)) {
        MSigDB.mm <- msigdbr::msigdbr(species = "Mus musculus",  # returned gene type : Entrez
                                      category = 'H')
        
        MSigDB.mm$gs_name <- sapply(MSigDB.mm$gs_name,
                                    function(x) strsplit(x, "HALLMARK_")[[1]][2])
        
    } else {
        MSigDB.mm <- msigdbr::msigdbr(species = "Mus musculus")  # returned gene type : Entrez
        MSigDB.mm <- MSigDB.mm %>% dplyr::filter(gs_name %in% target_gs)
    }
    
    MSigDB.mm.list <- split(MSigDB.mm$entrez_gene,  # each geneset is splitted by it's geneset name (returned as list)
                            MSigDB.mm$gs_name)
    
    gsva.data <- X_vst %>%
        SummarizedExperiment::assay() %>%  # returns matrix
        as.data.frame() %>%  # make dataframe before using `tibble` package
        tibble::rownames_to_column(var = "EnsemblID")  # default Gene ID of Salmon output : Ensembl ID
    
    gsva.data.IDmapped <- data.frame(
        "Entrez" = AnnotationDbi::mapIds(x = org.Mm.eg.db,
                                         keys = gsva.data$EnsemblID,
                                         keytype = "ENSEMBL",  # list of available keytypes can be found by `columns(org.Mm.eg.db)`
                                         column = "ENTREZID",  # list of available columns can be found by `columns(org.Mm.eg.db)`
                                         multiVals = "first")
    ) %>%
        dplyr::filter(!is.na(Entrez)) %>%  # remove unmapped genes
        tibble::rownames_to_column(var = "Ensembl") %>%  # rownames : Ensembl ID (key ID)
        dplyr::inner_join(gsva.data, by = c("Ensembl" = "EnsemblID"))  # merge with original count data
    
    gsva.data.geneMeans <- rowMeans(gsva.data.IDmapped %>% dplyr::select(-Ensembl, -Entrez))  # multi-mapping handling
    
    gsva.data.IDmapped <- gsva.data.IDmapped %>%  # multi-mapping handling
        dplyr::mutate(GeneMeans = gsva.data.geneMeans) %>%
        dplyr::select(Ensembl, Entrez, GeneMeans, dplyr::everything()) %>%
        dplyr::arrange(desc(GeneMeans)) %>%
        dplyr::distinct(Entrez, .keep_all = T) %>%  # if Entrez are same, retain whose value is the biggest
        dplyr::select(-Ensembl, -GeneMeans) %>%  # no need Ensembl ID and GeneMeans anymore
        tibble::column_to_rownames(var = "Entrez") %>%
        as.matrix()  # GSVA needs matrix
    
    gsva.param <- GSVA::gsvaParam(exprData = gsva.data.IDmapped,
                                  geneSets = MSigDB.mm.list,
                                  minSize = 5,
                                  maxSize = Inf,
                                  kcdf = "Gaussian",  # suit for vst transformed data
                                  tau = 1,  # calculate Gaussian-distributed scores
                                  maxDiff = T)
    
    gsva.result <- GSVA::gsva(param = gsva.param,
                              verbose = T)
    
    return(gsva.result)
}



#=======================================================================#
#   DDSAddSymbols : Add various types of gene symbols to DESeq result   #
#=======================================================================#
DDSAddSymbols <- function(dds_res_df){
    dds_df <- dds_res_df %>%
        as.data.frame() %>%  # in case when dds_res_df is not a dataframe
        tibble::rownames_to_column(var = "ENSG") %>%
        dplyr::inner_join(biomart.cache, by = c("ENSG" = "ensembl_gene_id")) %>%
        dplyr::select(ENSG, dplyr::any_of(biomart.cache %>% colnames()), dplyr::everything()) %>%  # https://stackoverflow.com/questions/22028937/
        dplyr::select(-ensembl_transcript_id) %>%
        tidyr::drop_na(ENSG, external_gene_name, entrezgene_id) %>%
        dplyr::distinct(ENSG, external_gene_name, entrezgene_id, .keep_all = T)
    
    return(dds_df)
}



#=======================================================================#
#   DDSCleanGenes : Remove unnamed/uninformed genes from DESeq result   #
#=======================================================================#
DDSCleanGenes <- function(dds_res_df) {
    
}



#=====================================#
#   Volcanoplot : draw Volcano plot   #
#=====================================#
Volcanoplot <- function(res_df, 
                        x = "log2FoldChange",
                        y = "padj",
                        gene_column = "external_gene_name",
                        target_genes = NULL,
                        title = "", 
                        subtitle = "",
                        caption = "",
                        legend_position = "none",
                        point_size = 3, 
                        label_size = 6,
                        pval_cutoff = 0.05,
                        lfc_cutoff = 0.58,
                        xlim = c(min(res_df[[x]], na.rm = TRUE) - 1.5,
                                 max(res_df[[x]], na.rm = TRUE) + 1.5),
                        ylim = c(0,
                                 max(-log10(res_df[[y]]), na.rm = TRUE) + 5)
) {
    # res_df shouldn't be filtered by pvalue or fold change
    
    EnhancedVolcano(res_df,
                    lab = res_df[[gene_column]],
                    selectLab = target_genes,
                    x = x,
                    y = y,
                    title = title,
                    subtitle = subtitle,
                    caption = caption,
                    legendPosition = legend_position,
                    pCutoff = pval_cutoff,
                    FCcutoff = lfc_cutoff,
                    pointSize = point_size,
                    labSize = label_size,
                    xlim = xlim,
                    ylim = ylim,
                    drawConnectors = T)
}



#===============================================================================#
#   enr_entrez2name : change entrez gene names to gene symbol in GSEA results   #
#===============================================================================#
enr_entrez2name <- function(gene_chunk) {
    if (stringr::str_length(gene_chunk) == 1)
        return(gene_chunk)
    
    gene_vec <- gene_chunk %>% 
        strsplit('/') %>%
        unlist() %>%
        as.integer()
    
    mapped_vec <- biomart.cache %>%
        dplyr::filter(entrezgene_id %in% gene_vec) %>%
        dplyr::pull(external_gene_name) %>%
        unique() %>%
        paste(collapse = ' / ')
    
    return(mapped_vec)
}



#=========================================================#
#   DDSEnrich : run Enrichment tests using DESeq result   #
#=========================================================#
DDSEnrich <- function(dds_res_df, seed = 1123) {
    # colnames(res_df) should contain "log2FoldChange" and "entrezgene_id"
    # colnames(background_gene_df) should contain "entrezgene_id"
    
    background_gene_df <- dds_res_df %>% 
        dplyr::select(entrezgene_id) %>%
        dplyr::distinct(entrezgene_id)
    
    sig_res <- dds_res_df %>% 
        dplyr::filter(abs(log2FoldChange) > 0.58 & padj < 0.05)
    
    ego_up <- clusterProfiler::enrichGO(
        gene = sig_res %>% dplyr::filter(log2FoldChange > 0) %>% pull(entrezgene_id),
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",  # from keytypes(org.Mm.eg.db)
        ont = "BP",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe = background_gene_df$entrezgene_id %>% as.character(),
        qvalueCutoff = 0.05,
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE,
        pool = FALSE
    ) %>% as.data.frame() %>% dplyr::filter(Count > 3)
    
    ego_down <- clusterProfiler::enrichGO(
        gene = sig_res %>% dplyr::filter(log2FoldChange < 0) %>% pull(entrezgene_id),
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",  # from keytypes(org.Mm.eg.db)
        ont = "BP",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe = background_gene_df$entrezgene_id %>% as.character(),
        qvalueCutoff = 0.05,
        minGSSize = 5,
        maxGSSize = 500,
        readable = FALSE,
        pool = FALSE
    ) %>% as.data.frame() %>% dplyr::filter(Count > 3)
    
    ekegg_up <- clusterProfiler::enrichKEGG(
        gene = sig_res %>% dplyr::filter(log2FoldChange > 0) %>% pull(entrezgene_id),
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe = background_gene_df$entrezgene_id %>% as.character(),
        minGSSize = 5,
        maxGSSize = 500,
        qvalueCutoff = 0.05,
        use_internal_data = FALSE
    ) %>% as.data.frame() %>% dplyr::filter(Count > 3)
    
    ekegg_down <- clusterProfiler::enrichKEGG(
        gene = sig_res %>% dplyr::filter(log2FoldChange < 0) %>% pull(entrezgene_id),
        organism = "mmu",
        keyType = "kegg",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        universe = background_gene_df$entrezgene_id %>% as.character(),
        minGSSize = 5,
        maxGSSize = 500,
        qvalueCutoff = 0.05,
        use_internal_data = FALSE
    ) %>% as.data.frame() %>% dplyr::filter(Count > 3)
    
    # https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
    FC_ordered_genes <- dds_res_df %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::pull(entrezgene_id)
    FC_ordered_vec <- dds_res_df %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::pull(log2FoldChange)
    names(FC_ordered_vec) <- FC_ordered_genes
    FC_ordered_vec <- na.omit(FC_ordered_vec)
    
    MSigDB.hallmark.mm <- msigdbr::msigdbr(species = "Mus musculus",  # returned gene type : Entrez
                                           category = 'H')
    
    MSigDB.hallmark.mm$gs_name <- sapply(MSigDB.hallmark.mm$gs_name,
                                         function(x) strsplit(x, "HALLMARK_")[[1]][2])
    
    MSigDB.hallmark.mm <- MSigDB.hallmark.mm %>% dplyr::select(gs_name, entrez_gene)
    
    gsea <- clusterProfiler::GSEA(
        geneList = FC_ordered_vec,
        TERM2GENE = MSigDB.hallmark.mm,
        minGSSize = 5,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        seed = seed
    ) %>% 
        as.data.frame() %>% 
        dplyr::select(-Description) %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        dplyr::mutate(core_enrichment = sapply(core_enrichment, enr_entrez2name))
    
    # only terms with p.adjust less than 0.05 will be returned
    gsego <- clusterProfiler::gseGO(
        geneList = FC_ordered_vec,
        ont = "BP",
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        exponent = 1,
        minGSSize = 5,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = seed,
        by = "fgsea"
    ) %>% as.data.frame() %>% dplyr::arrange(desc(NES))
    
    # only terms with p.adjust less than 0.05 will be returned
    gsekegg <- clusterProfiler::gseKEGG(
        geneList = FC_ordered_vec,
        organism = "mmu",
        keyType = "kegg",
        exponent = 1,
        minGSSize = 5,
        maxGSSize = 500,
        eps = 1e-10,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = TRUE,
        use_internal_data = FALSE,
        seed = seed,
        by = "fgsea"
    ) %>% as.data.frame() %>% dplyr::arrange(desc(NES))
    
    res <- list(ego_up, ego_down, ekegg_up, ekegg_down, gsea, gsego, gsekegg)
    names(res) <- c("ego_up", "ego_down", "ekegg_up", "ekegg_down", "gsea", "gsego", "gsekegg")
    
    return(res)
}