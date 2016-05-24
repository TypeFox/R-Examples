# Delimiters used in PED files.
delims <- "[ \t]"

#' Creates a matrix wrapper around binary PED files.
#'
#' \code{BEDMatrix} is an S3 class that behaves similarly to a regular
#' \code{matrix} by implementing key methods such as \code{[}, \code{dim}, and
#' \code{dimnames}. Subsets are extracted directly and on-demand from the binary
#' PED file without loading the entire file into memory through memory mapping.
#'
#' A \code{BEDMatrix} instance can be created by providing the path to the BED
#' file (with or without extension) as \code{path}, the number of individuals
#' as \code{n}, and the number of markers as \code{p}. If a FAM file (which
#' corresponds to the first six columns of a PED file) of the same name and in
#' the same directory as the BED file exists, it is optional to provide
#' \code{n} and the number of individuals as well as the rownames of the
#' \code{BEDMatrix} will be detected automatically. If a BIM file (which
#' corresponds to the MAP file that accompanies a PED file) of the same name
#' and in the same directory as the BED file exists, it is optional to provide
#' \code{p} and the number of markers as well as the colnames of the
#' \code{BEDMatrix} will be detected automatically.  For very large BED file it
#' is advised to provide \code{n} and \code{p} manually to speed up object
#' creation. In that case \code{rownames} and \code{colnames} will be set to
#' \code{NULL} and have to be specified manually.
#'
#' A BED file can be created from a PED file with
#' \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK} using \code{plink
#' --file myfile --make-bed}. BED files are storage and query efficient, and can
#' be transformed back into the original PED file with PLINK using \code{plink
#' --bfile myfile --recode}.
#'
#' Internally, \code{BEDMatrix} inherits from \code{list} and
#' exposes a few attributes that should not be relied upon in actual code:
#' \code{path}, \code{dims}, \code{dnames}, and \code{_instance}. \code{path}
#' stores the path to the BED file. \code{dims} and \code{dnames} contain the
#' dimensions and dimnames of the BEDMatrix object. \code{_instance} points to
#' the underlying \code{Rcpp} module. The \code{Rcpp} module exposes an S4 class
#' called \code{BEDMatrix_} that memory maps the BED file via
#' \code{Boost.Interprocess} of the \code{BH} package.
#'
#' @param path Path to the binary PED file, with or without extension.
#' @param n The number of individuals. Optional if FAM file of same name as BED
#'   file exists. If provided, \code{rownames} will be set to \code{NULL} and
#'   have to be provided manually.
#' @param p The number of markers. Optional if MAP file of same name as BED file
#'   exists. If provided, \code{colnames} will be set to \code{NULL} and have
#'   to be provided manually.
#' @export
BEDMatrix <- function(path, n = NULL, p = NULL) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        # Try to add extension (common in PLINK).
        path <- paste0(path, ".bed")
        if (!file.exists(path)) {
            stop("File not found.")
        }
    }
    dir <- substr(path, 1, nchar(path) - 4)
    if (is.null(n)) {
        # Check if FAM file exists.
        if (!file.exists(paste0(dir, ".fam"))) {
            stop("FAM file of same name not found. Provide number of individuals (n).")
        } else {
            message("Extracting number of individuals and rownames from FAM file...")
            fam <- readLines(paste0(dir, ".fam"))
            # Determine n.
            n <- length(fam)
            # Determine rownames.
            rownames <- sapply(strsplit(fam, delims), function(line) {
                # Concatenate family ID and subject ID
                return(paste0(line[1], "_", line[2]))
            })
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if BIM file exists.
        if (!file.exists(paste0(dir, ".bim"))) {
            stop("BIM file of same name not found. Provide number of markers (p).")
        } else {
            message("Extracting number of markers and colnames from BIM file...")
            bim <- readLines(paste0(dir, ".bim"))
            # Determine p.
            p <- length(bim)
            # Determine colnames.
            colnames <- sapply(strsplit(bim, delims), function(line) {
                # Concatenate SNP name and minor allele (like --recodeA)
                return(paste0(line[2], "_", line[5]))
            })
        }
    } else {
        p <- as.integer(p)
        colnames <- NULL
    }
    # Create Rcpp object
    rcpp_obj <- new(BEDMatrix_, path, n, p)
    # Wrap object in S3 class
    s3_obj <- list()
    class(s3_obj) <- "BEDMatrix"
    attr(s3_obj, "_instance") <- rcpp_obj
    attr(s3_obj, "dnames") <- list(rownames, colnames)
    attr(s3_obj, "dims") <- c(rcpp_obj$n, rcpp_obj$p)
    attr(s3_obj, "path") <- path
    return(s3_obj)
}

#' @export
print.BEDMatrix <- function(x, ...) {
    dims <- dim(x)
    n <- dims[1]
    p <- dims[2]
    cat(paste(n, "x", p, "BEDMatrix"), "\n")
}

#' @export
`[.BEDMatrix` <- function(x, i, j, drop = TRUE) {
    rcpp_obj <- attr(x, "_instance")
    dims <- dim(x)
    n <- dims[1]
    p <- dims[2]
    if (nargs() > 2) {
        # Case [i, j]
        if (missing(i)) {
            i <- 1:n
        } else if (class(i) == "logical") {
            i <- which(rep_len(i, n))
        } else if (class(i) == "character") {
            i <- sapply(i, function(name) {
                which(rownames(x) == name)
            }, USE.NAMES = FALSE)
        }
        if (missing(j)) {
            j <- 1:p
        } else if (class(j) == "logical") {
            j <- which(rep_len(j, p))
        } else if (class(j) == "character") {
            j <- sapply(j, function(name) {
                which(colnames(x) == name)
            }, USE.NAMES = FALSE)
        }
        subset <- rcpp_obj$matrixSubset(x, i, j)
        # Let R handle drop behavior.
        if (drop == TRUE && (nrow(subset) == 1 || ncol(subset) == 1)) {
            subset <- subset[, ]
        }
    } else {
        if (missing(i)) {
            # Case []
            i <- 1:n
            j <- 1:p
            subset <- rcpp_obj$matrixSubset(x, i, j)
        } else {
            # Case [i]
            if (class(i) == "matrix") {
                i <- as.vector(i)
                if (class(i) == "logical") {
                  i <- which(rep_len(i, n * p))
                  # matrix treats NAs as TRUE
                  i <- sort(c(i, which(is.na(x[]))))
                }
            } else {
                if (class(i) == "logical") {
                  i <- which(rep_len(i, n * p))
                }
            }
            subset <- rcpp_obj$vectorSubset(x, i)
        }
    }
    return(subset)
}

#' @export
dim.BEDMatrix <- function(x) {
    attr(x, "dims")
}

#' @export
dimnames.BEDMatrix <- function(x) {
    attr(x, "dnames")
}

#' @export
`dimnames<-.BEDMatrix` <- function(x, value) {
    d <- dim(x)
    v1 <- value[[1]]
    v2 <- value[[2]]
    if (!is.list(value) || length(value) != 2 || !(is.null(v1) || length(v1) == d[1]) || !(is.null(v2) || length(v2) == d[2])) {
        stop("invalid dimnames")
    }
    attr(x, "dnames") <- lapply(value, function(v) {
        if (!is.null(v)) {
            as.character(v)
        }
    })
    return(x)
}

#' @export
length.BEDMatrix <- function(x) {
    prod(dim(x))
}

#' @export
is.matrix.BEDMatrix <- function(x) {
    TRUE
}

#' @export
as.matrix.BEDMatrix <- function(x, ...) {
    x[, , drop = FALSE]
}
