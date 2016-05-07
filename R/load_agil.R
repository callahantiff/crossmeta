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

load_agil <- function (gse_names, data_dir) {

    esets <- list()
    for (gse_name in gse_names) {

        gse_dir <- file.path(data_dir, gse_name)
        save_name <- paste(gse_name, "eset.rds", sep="_")
        eset_path <- list.files(gse_dir, save_name, full.names=TRUE)

        #check if saved copy
        if (length(eset_path) != 0) {
            eset <- readRDS(eset_path)

        } else {
            #get GSEMatrix (for pheno data)
            eset <- GEOquery::getGEO(gse_name, destdir=gse_dir, GSEMatrix=TRUE)[[1]]

            #load non-normalized txt files and normalize
            data_paths <- list.files(gse_dir, pattern="GSM.*txt",
                                     full.names=TRUE, ignore.case=TRUE)
            data <- limma::read.maimages(data_paths,
                                         source="agilent", green.only=TRUE)
            data <- limma::neqc(data, status=data$genes$ControlType,
                                negctrl=-1, regular=0)

            #fix up sample/feature names
            colnames(data) <- stringr::str_match(colnames(data),
                                                 ".*(GSM\\d+).*")[, 2]
            row.names(data$E) <- data$genes$ProbeName
            data$genes <- data$genes[!duplicated(data$genes$ProbeName),]
            row.names(data$genes) <- data$genes$ProbeName

            #transfer to eset
            exprs(eset) <- data$E
            fData(eset) <- data$genes[featureNames(eset), ]

            #add SYMBOL annotation
            gpl_name <- annotation(eset)
            eset <- symbol_annot(eset, gpl_name)

            #save to disc
            saveRDS(eset, file.path(gse_dir, save_name))
        }
        esets[[gse_name]] <- eset
    }
    return(esets)
}