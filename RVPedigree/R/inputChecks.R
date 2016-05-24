## This file contains checks of input function parameters that are
## used by several other functions.

check_method <- function(method) {
    if (! method %in% c("ASKAT", "NASKAT", "VCC1", "VCC2", "VCC3")) {
        stop("Unkown method specified; the method parameter ",
             "should be one of: 'ASKAT', 'NASKAT', 'VCC1', ",
             "'VCC2', or 'VCC3'" )
    }
}


check_pheno <- function(y) {
    if (is.null(y)) {
        stop("Argument y is empty, you did not specify a phenotype vector")
    }

    if (class(y) == "matrix") {
        if(ncol(y) == 1){
            y <- as.vector(y)
        }
    }

    if (!class(y) == "numeric") {
        stop("Argument y should be a numeric vector")
    }

    return(y)
}


check_covariates <- function(X, y) {
    if (!is.null(X)) {
        if (!class(X) == "matrix") {
            stop("Argument X is not a matrix")
        }
    }

    if (nrow(X) != length(y)) {
        msg <- paste0("The number of rows of the covariate matrix (",
                      nrow(X),
                      ") is not equal to the number of elements in the",
                      " phenotype vector (",
                      length(y),
                      ")")
        stop(msg)
    }
}


check_relmatrix <- function(Phi) {
    if (!class(Phi) == "matrix") {
        stop("Argument Phi is not a matrix")
    }
}


check_files <- function(filename, type) {
    if (is.null(type)) {
        stop("Argument type was not specified")
    }

    if (is.null(filename)) {
        stop("No file name for genotype data was specified")
    }

    if (!file.exists(filename)) {
        stop(paste("File", filename, "does not exist"))
    }

    fileExtension <- substr(filename,
                            nchar(filename) - 3,
                            nchar(filename))
    if (fileExtension == ".ped" & type == "bed") {
        warning("The file name for the genotype data ends in 'ped', ",
                "but you specified type='bed'; reading genotype data",
                " will probably fail")
    }
    if (fileExtension == ".bed" & type == "ped") {
        warning("The file name for the genotype data ends in 'bed', ",
                "but you specified type='ped'; reading genotype data",
                " will probably fail")
    }
}


check_positions <- function(startpos, endpos) {
    if (startpos > endpos) {
        stop("The start bp position is larger than the end bp position")
    }

    if (startpos == endpos) {
        stop("The start bp position is equal to the end bp position")
    }

    if (startpos < 0 | endpos < 0) {
        stop("The startpos and endpos parameters should be >= 0")
    }
}


check_weights <- function(weights) {
    if (!is.null(weights)) {
        if(!is.vector(weights) | class(weights) != "numeric") {
            stop("The 'weights' parameter should be a numeric vector")
        }
    }
}


check_regions <- function(regions) {
    if (is.null(regions)) {
        stop("No data frame with genomic regions was specified.")
    }

    if (class(regions) != "data.frame") {
        stop("The 'regions' parameter should be a data frame")
    }

    if (ncol(regions) < 4) {
        msg <- paste("The 'regions' parameter should have at least 4",
                     " columns: 'Name', 'Chr', 'StartPos' and 'EndPos'")
        stop(msg)
    }

    columnNames <- colnames(regions)
    if (! "Name" %in% columnNames) {
        stop("The 'regions' data frame doesn't contain a column named 'Name'")
    }

    if (! "Chr" %in% columnNames) {
        stop("The 'regions' data frame doesn't contain a column named 'Chr'")
    }

    if (! "StartPos" %in% columnNames) {
        stop("The 'regions' data frame doesn't contain a column named 'StartPos'")
    }

    if (! "EndPos" %in% columnNames) {
        stop("The 'regions' data frame doesn't contain a column named 'EndPos'")
    }
}


check_Ncores <- function(Ncores) {
        if (Ncores < 1) {
        stop("The number of CPU cores must be >= 1")
    }
}


check_VCC3afterVCC1 <- function(VCC3afterVCC1, method) {
    if (VCC3afterVCC1) {
        if (method != "VCC1") {
            warning("'VCC3afterVCC1' parameter set, but method",
                    " is not 'VCC1'; ignoring this parameter")
        }
    }
}
