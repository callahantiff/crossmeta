#########################################################################################
### load_agil - modified code from crossmeta that loads and processes agilent data
### version 1.0.0
### date: 03.23.18
#########################################################################################

source("~/Dropbox/GraduateSchool/PhD/LabWork/MetaOmic/transcriptomic/R/crossmeta/R/load_illum.R")

# ----
# Load and pre-process raw Agilent files.
#
# Load raw txt files previously downloaded with \code{get_raw_agil}. Used by
# \code{load_raw}.
#
# Data is normalized and SYMBOL and PROBE annotation are added to fData slot.
#
# @param gse_names Character vector of Agilent GSE names.
# @param data_dir String specifying directory with GSE folders.
#
# @seealso \code{\link{get_raw}} to obtain raw data.
# @return List of annotated esets.

load_agil <- function (gse_names, data_dir, gpl_dir, ensql) {
    require(Biobase, quietly = TRUE)
  
    esets  <- list()
    errors <- c()
    
    for (gse_name in gse_names) {
        cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Processing AGILENT:", gse_name, strrep("#", 5), "\n")
        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")
        
        ## check for raw data files
        # ensure that we have non-normalized data
        file_paths <- list.files(gse_dir, "*[^matrix].txt", full.names = TRUE, ignore.case = TRUE)
        if (length(file_paths) < 2) stop(cat(gse_name, ': no raw data provided.'))

        # get GSEMatrix (for pheno data)
        eset <- NULL
        while (is.null(eset)) {
          eset <- try(suppressWarnings(getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE)))
        }

        # check if have GPL
        cat(strrep("#", 5), "Downloading Annotation Files", strrep("#", 5), "\n")
        gpl_names <- paste0(sapply(eset, Biobase::annotation), '.soft', collapse = "|")
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
            list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]})

        # copy over GPL
        if (length(gpl_paths) > 0) file.copy(gpl_paths, gse_dir)

        # will use local GPL or download if couldn't copy
        eset <- NULL
        while (is.null(eset)) {
            eset <- try(suppressWarnings(getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE)))
        }

        # name esets
        if (length(eset) > 1) {
            names(eset) <- paste(gse_name, sapply(eset, annotation), sep='.')
        } else {
            names(eset) <- gse_name
        }
        
        # determine is single or double channel
        pdat_names <- sapply(eset, function(x) colnames(Biobase::pData(x)))
        ch = ifelse((any(grepl('ch2', pdat_names))), 'two', 'one')
        cat(gse_name, ": is a", ch, "channel array", "\n")

        # load eset for each platform in GSE
        cat(strrep("#", 5), "Loading Expression Data", strrep("#", 5), "\n")
        eset <- lapply(eset, function(eset.gpl) {
            tryCatch(load_agil_plat(eset.gpl, ch, gse_dir, gse_name, ensql),
                     error = function(e) NA)
        })

        # save to disc
        if (!all(is.na(eset)))
            saveRDS(eset[!is.na(eset)], file.path(gse_dir, save_name))
        if (anyNA(eset))
            errors <- c(errors, names(eset[is.na(eset)]))
        if (!all(is.na(eset)))
            esets[[gse_name]] <- eset[!is.na(eset)]
    }

    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (list(esets = esets, errors = errors))
}


# ----
# Helper utility for load_agil.
#
# This function loads the data, merges the expression data features with Elist of expression data,
# verifies log transformation, and performs RMA background correction (against control probes) and quantile normalization
#
# Used by load_agil to load an eset for each GPL platform in a GSE.
#
# @param eset Expression set obtained by load_agil call to getGEO.
# @param ch String specifying pthe number of channels.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_agil}}.
# @return Annotated eset with scan_date in pData slot.
load_agil_plat <- function (eset, ch, gse_dir, gse_name, ensql) {
    require(limma, quietly = TRUE)
    require(Biobase, quietly = TRUE)
    require(stringr, quietly = TRUE)
    require(arrayQualityMetrics, quietly = TRUE)

    # try(Biobase::fData(eset)[Biobase::fData(eset) == ""] <- NA)
    try(Biobase::fData(eset)[] <- lapply(Biobase::fData(eset), as.character))
  
    # define groups
    groups = grep("characteristics", names(Biobase::pData(eset)), value=TRUE)
    # study_name = ifelse(is.null(names(eset)), gse_name, names(eset))
    study_name = paste(gse_name, eset[1]@annotation, sep = "_")
  
    # # Quality control pre-processing
    # cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Quality Assessment: Pre-Processing", study_name, strrep("#", 5), "\n")
    # arrayQualityMetrics(expressionset = eset,
    #                     outdir = paste(gse_dir, "/QualityControl/", study_name, "_QualityControl_Report_Pre", sep = ""),
    #                     reporttitle = as.character(paste(study_name, "Pre-Processing Quality Control Report")),
    #                     intgroup = c(groups, "geo_accession"),
    #                     do.logtransform = TRUE,
    #                     force = TRUE)

    # get paths to raw files for samples in eset
    pattern <- paste(sampleNames(eset), ".*txt", collapse = "|", sep = "")
    elist_paths <- list.files(gse_dir, pattern, full.names = TRUE, ignore.case = TRUE)

    # if multiple with same GSM, take first
    gsm_names  <- stringr::str_extract(elist_paths, "GSM[0-9]+")
    elist_paths <- elist_paths[!duplicated(gsm_names)]

    ## process differently depending on single vs. double channel array
    # single channel
    if (ch == 'one') {
      # load non-normalized txt files and normalize
      elist <- tryCatch(limma::read.maimages(elist_paths, source = "agilent", green.only = TRUE),
                        error = function(e) {
                          # determine source of error
                          output <- capture.output(tryCatch(
                            limma::read.maimages(elist_paths, source = "agilent", green.only = TRUE),
                            error = function(e) NULL))
                          
                          # retry with error excluded
                          exclude     <- which(elist_paths == gsub('^Read| ', '', output[length(output)-1])) + 1
                          elist_paths <- elist_paths[-exclude]
                          limma::read.maimages(elist_paths, source = "agilent", green.only = TRUE)
                        })
      
      # rma background correction and quantile normalize performed at the same time data is log2 transformed
      elist <- limma::neqc(elist, status = elist$genes$ControlType, negctrl = -1, regular = 0)
      
      # fix up sample names
      colnames(elist) <- stringr::str_match(colnames(elist), ".*(GSM\\d+).*")[, 2]
      eset <- eset[, colnames(elist)]
      
      # merge elist and eset feature data
      elist <- merge_elist(eset, elist)
      row.names(elist$E) <- row.names(elist$genes) <- make.unique(elist$genes$ProbeName)
      
      # transfer to eset
      eset <- ExpressionSet(elist$E,
                            phenoData = Biobase::phenoData(eset),
                            featureData = as(elist$genes, 'AnnotatedDataFrame'),
                            annotation = Biobase::annotation(eset))
    }
    
    # double channels
    # (http://matticklab.com/index.php?title=Two_channel_analysis_of_Agilent_microarray_data_with_Limma)
    if (ch == 'two') {
      # load non-normalized txt files and normalize
      elist <- tryCatch(limma::read.maimages(elist_paths, source = "agilent"),
                        error = function(e) {
                          # determine source of error
                          output <- capture.output(tryCatch(
                            limma::read.maimages(elist_paths, source = "agilent"),
                            error = function(e) NULL))
                          
                          # retry with error excluded
                          exclude     <- which(elist_paths == gsub('^Read| ', '', output[length(output)-1])) + 1
                          elist_paths <- elist_paths[-exclude]
                          limma::read.maimages(elist_paths, source = "agilent")
                        })
      
      # rma background correction and quantile normalize performed at the same time data is log2 transformed
      elist <- limma::backgroundCorrect(elist, method = "normexp", normexp.method = "rma")
      elist <- limma::normalizeWithinArrays(elist, method = "loess")
      elist <- limma::normalizeBetweenArrays(elist, method = "quantile")
      
      # fix up sample names
      colnames(elist) <- stringr::str_match(colnames(elist), ".*(GSM\\d+).*")[, 2]
      eset <- eset[, colnames(elist)]

      # merge elist and eset feature data
      elist <- merge_elist(eset, elist)
      rownames(elist$M) <- rownames(elist$genes) <- make.unique(elist$genes$ProbeName)
      
      # transfer to eset - making the assumption that all of the two-color microarray experiments will have used a common reference for all microarrays
      eset <- ExpressionSet(elist$M,
                            phenoData = Biobase::phenoData(eset),
                            featureData = as(elist$genes, 'AnnotatedDataFrame'),
                            annotation = Biobase::annotation(eset))
    }
    
    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name, ensql)

    # Quality control post-processing
    cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Quality Assessment: Post-Processing", study_name, strrep("#", 5), "\n")
    arrayQualityMetrics(expressionset = eset,
                        outdir = paste(gse_dir, "/QualityControl/", study_name, "_QualityControl_Report_Post", sep = ""),
                        reporttitle = as.character(paste(study_name, "Post-Processing Quality Control Report")),
                        intgroup = c(groups, "geo_accession"),
                        do.logtransform = FALSE,
                        force = TRUE)
  
    return(eset)
}

