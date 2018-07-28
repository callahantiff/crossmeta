#########################################################################################
### load_illum - modified code from crossmeta that loads and processes illumina data
### version 1.0.0
### date: 03.20.18
#########################################################################################

# load needed scripts, can remove this again once we have packaged the code
source("~/Dropbox/GraduateSchool/PhD/LabWork/MetaOmic/transcriptomic/R/crossmeta/R/illum_headers.R")


# ----
# Load and pre-process raw Illum files.
#
# Load raw txt files previously downloaded with \code{get_raw} and checked
# for format with \code{open_raw_illum}. Used by \code{load_raw}.
#
# Data is normalized, SYMBOL and PROBE annotation are added to fData slot, and
# detection p-values are added to pvals slot.
#
# @param gse_names Character vector of Illumina GSE names.
# @param data_dir String specifying directory with GSE folders.
#
# @seealso \code{\link{get_raw}} to obtain raw Illumina data.
#   \code{\link{open_raw_illum}} to ensure their correctness.

# @return List of annotated esets.

load_illum <- function (gse_names, data_dir, gpl_dir, ensql) {
  require(Biobase, quietly = TRUE)
    esets  <- list()
    errors <- c()
    gse_meta <- data.frame(study = character(),
                           platform = character(),
                           channels = character(),
                           data_type = character(),
                           sample_number = numeric(),
                           samples = character(),
                           annotation_package = character(),
                           num_probe = numeric(),
                           num_annotated_probes_all = numeric(),
                           num_annotated_probes_unique = numeric(),
                           stringsAsFactors = FALSE)
    
    for (gse_name in gse_names) {
      cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Processing ILLUMINA:", gse_name, strrep("#", 5), "\n")
        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep = "_")

        # get GSEMatrix (for pheno dat)
        eset <- NULL
        while (is.null(eset)) {
          eset <- try(suppressWarnings(getGEO(gse_name, destdir = gse_dir, GSEMatrix = TRUE, getGPL = FALSE)))
        }
        
        # check if have GPL
        cat(strrep("#", 5), "Downloading Annotation Files", strrep("#", 5), "\n")
        gpl_names <- paste0(sapply(eset, Biobase::annotation), '.soft', collapse = "|")
        gpl_paths <- sapply(gpl_names, function(gpl_name) {
          list.files(gpl_dir, gpl_name, full.names = TRUE, recursive = TRUE, include.dirs = TRUE)[1]
        })
        
        # copy over GPL
        if (length(gpl_paths) > 0) file.copy(gpl_paths, gse_dir)
        
        # if (length(eset) > 1) {
        #     warning("Multi-platform Illumina GSEs not supported. ", gse_name)
        #     errors <- c(errors, gse_name)
        #     next
        # }
        
        # name esets
        if (length(eset) > 1) {
          names(eset) <- paste(gse_name, sapply(eset, Biobase::annotation), sep = '.')
        } else {
          names(eset) <- gse_name
        }
        
        # load eset for each platform in GSE
        eset <- lapply(eset, function(set) {
          cat(strrep("#", 5), "Loading Expression Data", strrep("#", 5), "\n")
          tryCatch(
            load_illum_plat(set, gse_dir, gse_name, ensql),
            error = function(e)
              NULL)
        })
        
        # save to disc
        if (anyNA(eset))
          errors <- c(errors, names(eset[is.na(eset)]))
        
        if (!all(is.na(eset))) {
          # esets[[gse_name]] <- eset[!is.na(eset)]
          saveRDS(eset[!is.na(eset)], file.path(gse_dir, save_name))
          
          if (length(eset) > 1) {
            for(i in 1:length(eset)) {
              esets[[names(eset[i])]] <- eset[i][!is.na(eset[i])]
            }
          }
          if (length(eset) == 1) {
            esets[[gse_name]] <- eset[!is.na(eset)]
          }
        }
    }
    
    # compile metadata
    cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Compiling Study Metadata", strrep("#", 5), "\n")
    for(i in 1:length(esets)){
      if (length(esets[[i]]) > 0) {
        res = esets[[i]]
        gse_meta[i, "study"] = names(res)
        gse_meta[i, "platform"] = gpl_bioc[which(rownames(gpl_bioc) == Biobase::annotation(res[[1]])),]$title
        gse_meta[i, "channels"] = "one"
        gse_meta[i, "data_type"] = ifelse(length(grep("non.*norm.*txt$|raw.*txt$|nonorm.*txt*", 
                                                      phenoData(res[[1]])$supplementary_file)) == length(sampleNames(res)), 
                                          "Raw Microarray Data",
                                          ifelse(length(grep("*non-normalized", phenoData(res[[1]])$supplementary_file)) > 0,
                                                 "Non-Normalized Gene Expression Matrix",
                                                 "Normalized Gene Expression Matrix"))
        gse_meta[i, "sample_number"] = length(sampleNames(res))
        gse_meta[i, "samples"] = paste(unlist(sampleNames(res)), collapse=',')
        gse_meta[i, "annotation_package"] = ifelse(gpl_bioc[which(rownames(gpl_bioc) == Biobase::annotation(res[[1]])),]$bioc_package == "",
                                                   "AnnotationDbi::species(biocpack)",
                                                   gpl_bioc[which(rownames(gpl_bioc) == Biobase::annotation(res[[1]])),]$bioc_package)
        gse_meta[i, "num_probe"] = length(fData(res[[1]])$PROBE)
        gse_meta[i, "num_annotated_probes_all"] = length(na.omit(fData(res[[1]])$ENTREZID))
        gse_meta[i, "num_annotated_probes_unique"] = length(unique(na.omit(fData(res[[1]])$ENTREZID)))
      }
    }
    
    eset_names <- get_eset_names(esets, gse_names)
    esets <- unlist(esets)
    names(esets) <- eset_names
    return (list(esets = esets, errors = errors, metadata = gse_meta))
}

# Helper utility for load_illum.
#
# This function loads the data, merges the expression data features with Elist of expression data,
# verifies log transformation, and performs RMA background correction and quantile normalization
#
# Used by load_illum to load an eset.
#
# @param eset Expression set obtained by load_illum call to getGEO.
# @param gse_name String specifying GSE name.
# @param gse_dir String specifying path to GSE folder.
#
# @seealso \code{\link{load_illum}}.
# @return Annotated eset.

# ----
load_illum_plat <- function(eset, gse_dir, gse_name, ensql) {
    require(Biobase, quietly = TRUE)
    require(limma, quietly = TRUE)
    require(arrayQualityMetrics, quietly = TRUE)

    # try(Biobase::fData(eset)[Biobase::fData(eset) == ""] <- NA)
    try(Biobase::fData(eset)[] <- lapply(Biobase::fData(eset), as.character))
  
    # define groups - needed for checking post-processing array quality
    groups = grep("characteristics", names(Biobase::pData(eset)), value=TRUE)
    study_name = paste(gse_name, eset[1]@annotation, sep = "_")

    # fix header issues
    elist_paths <- list.files(gse_dir, pattern = "non.*norm.*txt$|raw.*txt$|nonorm.*txt$", full.names = TRUE, ignore.case = TRUE)
    
    # for multiplatform arrays
    if (length(elist_paths[!grepl('fixed[.]txt$', elist_paths)]) > 1) {
      elist_paths <- elist_paths[!grepl('fixed[.]txt$', elist_paths)]
      elist_paths <- elist_paths[grepl(as.character(Biobase::pData(eset)$description[[1]]), tolower(elist_paths))]
      annotation  <- fix_illum_headers(elist_paths, eset)
    }
    else {
      elist_paths <- elist_paths[!grepl('fixed[.]txt$', elist_paths)]
      annotation  <- fix_illum_headers(elist_paths, eset)
    }

    # load fixed elist paths
    elist_paths <- gsub(".txt", "_fixed.txt", elist_paths, fixed = TRUE)
    elist <- limma::read.ilmn(elist_paths, probeid = "ID_REF", annotation = annotation)
    
    # don't correct if already log transformed (already corrected?)
    logd <- max(elist$E, na.rm = TRUE) < 100
    if (!logd) {
        elist <- tryCatch (
            limma::neqc(elist),
            error = function(c) {
                # PMID:19068485 recommends mle and offset 50
                # elist <- limma::backgroundCorrect(elist, method = "normexp", normexp.method = "mle", offset = 50)
                elist <- limma::backgroundCorrect(elist, method = "normexp", normexp.method = "rma")
                return(limma::normalizeBetweenArrays(elist, method = "quantile"))
            })
    }

    # merge eset and elist fdata
    elist <- merge_elist(eset, elist)
    if ('ID_REF' %in% colnames(elist$genes)) {
        row.names(elist$E) <- row.names(elist$genes) <- make.unique(elist$genes$ID_REF)
    } else {
        row.names(elist$E) <- row.names(elist$genes) <- NULL
    }

    # determine best sample matches
    res   <- match_samples(eset, elist)
    elist <- elist[, res$elist_order]
    eset  <- eset[, res$eset_order]
    warn  <- res$warn

    # keep gse matrix and raw elist title
    Biobase::pData(eset)$title.gsemat   <- Biobase::pData(eset)$title
    Biobase::pData(eset)$title.raw      <- colnames(elist)

    if (warn) {
        # use raw elist titles to ensure correct contrasts
        Biobase::pData(eset)$title <- colnames(elist)

        # add illum colname to warn about pData
        Biobase::pData(eset)$illum <- NA
    }
    colnames(elist) <- Biobase::sampleNames(eset)

    # transfer elist to eset
    eset <- Biobase::ExpressionSet(elist$E,
                          phenoData = Biobase::phenoData(eset),
                          featureData = as(elist$genes, 'AnnotatedDataFrame'),
                          annotation = Biobase::annotation(eset))
    
    # add SYMBOL annotation
    eset <- symbol_annot(eset, gse_name, ensql)
    
    # quality control post-processing
    # cat("\n\n", strrep("#",100), "\n", strrep("#", 5), "Quality Assessment: Post-Processing", study_name, strrep("#", 5), "\n")
    # arrayQualityMetrics(expressionset = eset,
    #                     outdir = paste(gse_dir, "/QualityControl/", study_name, "_QualityControlReport", sep = ""),
    #                     reporttitle = as.character(paste(study_name, "Post-Processing Quality Control Report")),
    #                     intgroup = c(groups, "geo_accession"),
    #                     do.logtransform = FALSE,
    #                     force = TRUE)

    return(eset)
}

# -------------------

fuzzy_pmatch <- function(x, table) {
  # like base pmatch
  # partial match occurs if the whole of the element of x matches any part of the element of table
    x1 <- gsub( " |[[:punct:]]+", "", as.character(tolower(x)))
    table1 <- gsub( " |[[:punct:]]+", "", as.character(tolower(table)))

    # first look for perfect matches
    perfect <- match(x1, table1)

    # is every x has a perfect match in table, return
    if (sum(is.na(perfect)) == 0) {
      # table = x
      return(perfect)
    }

    # otherwise first grep
    tomatch <- x1[is.na(perfect)]
    gmatch  <- sapply(tomatch, function(val) {
        res <- grep(val, table1, fixed = TRUE)[1]
        if (!length(res)) return(NA_integer_)
        return(res)
    })

    # fill in grep result where NA in perfect
    perfect[is.na(perfect)] <- gmatch[is.na(perfect)]
    return(perfect)
}

# ----
match_samples <- function(eset, elist) {
  # function ensures that the same identifers for samples are used with eset and elist
    require(Biobase, quietly = TRUE)
    require(ccmap, quietly = TRUE)
  
    # determine if elist has fewer samples
    data_fewer <- ncol(elist) < ncol(eset)

    # check if colnames match
    if (data_fewer) {
        # check if all elist colnames in eset colnames
        if (all(colnames(elist) %in% colnames(eset))) {
            cat('Illumina samples matched by column names.\n')
            return(list(elist_order = colnames(elist), eset_order = colnames(elist), warn = FALSE))
        }
    } else {
        # check if all eset colnames in elist colnames
        if (all(colnames(eset) %in% colnames(elist))) {
            cat('Illumina samples matched by column names.\n')
            return(list(elist_order = colnames(eset), eset_order = colnames(eset), warn = FALSE))
        }
    }

    # check if eset pdata col matches elist colnames
    if (!is.null(colnames(elist))) {
        # matrix of positions of matches for elist colnames among those for each pdata column
        matches <- sapply(Biobase::pData(eset), function(col) {
            fuzzy_pmatch(colnames(elist), col)
        })

        # number of unique non NA matches for each pdata column
        nunique <- apply(matches, 2, function(match) length(unique(match[!is.na(match)])))

        # number of unique non NA matches should be the min of number of eset or pdata samples
        nmin <- min(ncol(eset), ncol(elist))
        if (any(nunique == nmin)) {
            cat('Illumina samples matched by pdata column.\n')

            # matches where satisfied
            bestcol <- names(which(nunique == nmin))[1]
            matches <- matches[, bestcol]

            # elist_order is positions where matches are not NA
            elist_order <- which(!is.na(matches))

            # eset_order is non NA matches
            eset_order <- matches[!is.na(matches)]

            return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))
        }
    }

    # ----
    # check if similarity offers unique match ----
    # make sure eset is log2 transformed
    logd <- max(Biobase::exprs(eset), na.rm = TRUE) < 1000
    if (!logd) {
      Biobase::exprs(eset) <- log2(Biobase::exprs(eset) + abs(min(Biobase::exprs(eset), na.rm = TRUE)) + 16)
    }

    # row names are the best match columns for elist and eset
    elist <- elist[!is.na(elist$genes[[1]]), ]
    row.names(elist) <- make.unique(elist$genes[[1]])
    # eset2 <- eset[!is.na(Biobase::fData(eset)), ]
    row.names(eset) <- make.unique(row.names(Biobase::fData(eset)))

    # only include rows without missing values
    eset   <- eset[complete.cases(Biobase::exprs(eset)), ]
    elist  <- elist[complete.cases(elist$E), ]
    qres   <- list()
    ngenes <- min(nrow(eset), nrow(elist))

    if (data_fewer) {
        # determine most similar eset sample for each sample in elist
        for (i in 1:ncol(elist)) {
            qsamp <- elist$E[, i]
            qres[[colnames(elist)[i]]] <- ccmap::query_drugs(qsamp, exprs(eset), sorted = FALSE, ngenes = ngenes)
        }

    } else {
        # determine most similar elist sample for each sample in eset
        for (i in 1:ncol(eset)) {
            qsamp <- Biobase::exprs(eset)[, i]
            qres[[colnames(eset)[i]]] <- ccmap::query_drugs(qsamp, elist$E, sorted = FALSE, ngenes = ngenes)
        }
    }

    # eset sample to most similar elist sample
    qres <- as.data.frame(qres)
    best <- sapply(qres, which.max)

    if (length(best) == length(unique(best))) {
        cat('Illumina samples matched by similarity.\n')

        if (data_fewer) {
            elist_order <- colnames(elist)
            eset_order <- best
        } else {
            elist_order <- best
            eset_order <- colnames(eset)
        }

        return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))

    } else {
        # look for misses in non-first query results
        dups   <- unique(best[duplicated(best)])
        misses <- setdiff(1:nrow(qres), unique(best))
        n <- nrow(qres)
        
        for (dup in dups) {
            # query results for duplicate
            i <- 1
            qres_dup  <- qres[, best == dup]

            while (dup %in% dups & i < n) {
                ibest_dup <- sapply(qres_dup, function(col) which(col == sort(col, partial=n-i)[n-i]))

                # for each miss
                for (miss in misses) {
                    # check if one ibest is miss
                    imiss <- ibest_dup == miss

                    if (sum(imiss) == 1){
                        # if so, replace best with ibest
                        ibest_repl <- ibest_dup[imiss]
                        best[names(ibest_repl)] <- ibest_repl

                        # also update duplicates and misses
                        dups   <- best[duplicated(best)]
                        misses <- setdiff(1:nrow(qres), unique(best))

                        # if no more misses, break
                        if (!length(misses)) {
                            break()
                        }
                    }
                }
                i <- i + 1
            }
        }

        if (!length(dups)) {
            cat('Illumina samples matched by similarity using non-first ranks.\n')
            if (data_fewer) {
                elist_order <- colnames(elist)
                eset_order <- best
            } else {
                elist_order <- best
                eset_order <- colnames(eset)
            }
            return(list(elist_order = elist_order, eset_order = eset_order, warn = FALSE))

        } else {
            cat('Illumina samples not matched.\n')
            return(list(elist_order = colnames(elist), eset_order = colnames(eset), warn = TRUE))
        }
    }
}

# ----
merge_elist <- function(eset, elist) {
  # function checks if there is feature data that comes with eset and if there is, it merges the
  # features by probe id with the elist data
  
    require(Biobase, quietly = TRUE)

    if (is.null(elist$genes)) stop('Raw elist lacks feature names.')

    # get eset and elist pdata columns
    esetcols  <- Biobase::fData(eset)
    elistcols <- elist$genes

    # first, check if the eset has feature data (not all do)
    if (ncol(esetcols) != 0) {
      # find eset pdata column that best matches elist features
      best  <- c(esetcol=NA, elistcol=NA)
      bestf <- 0
      
      for (i in seq_along(elistcols)) {
        elistcol <- elistcols[[i]]
        
        # get fraction of fdata column that has a match
        matches <- sapply(names(esetcols), function(esetcol) {
          sum(elistcol %in% esetcols[, esetcol]) / length(elistcol)
        })
        
        # update best
        if (max(matches) >= bestf) {
          bestf <- max(matches)
          best['elistcol'] <- names(elistcols)[i]
          best['esetcol']  <- names(matches[which.max(matches)])
        }
      }
      
      if (bestf > 0.3) {
        # merge eset and elist fdata columns
        esetcols  <- esetcols[!duplicated(esetcols[best['esetcol']]),, drop = FALSE]
        elistcols <- merge(elistcols, esetcols, all.x = TRUE, by.x = best['elistcol'], by.y = best['esetcol'], sort = FALSE)
        elistcols[elistcols == ""] <- NA
        elist$genes <- elistcols
        
        # add best info for illumina sample matching
        elist$elistcol <- best[['elistcol']]
        elist$esetcol  <- best[['esetcol']]
      }
    }
    
    elist$genes[] <- lapply(elist$genes, as.character)
    return(elist)
}

