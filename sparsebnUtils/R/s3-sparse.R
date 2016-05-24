#
#  s3-sparse.R
#  sparsebnUtils
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2016 Bryon Aragam. All rights reserved.
#

#------------------------------------------------------------------------------#
# sparse S3 Class for R
#------------------------------------------------------------------------------#

#
# An alternative data structure for storing sparse matrices in R using the (row, column, value)
#   format. Internally it is stored as a list with three components, each vectors, that contain
#   the rows / columns / values of the nonzero elements.
#
# Its main purpose is to serve as an intermediary between the standard R dense matrix class and the
#   internal SparseBlockMatrixR class. That is, to convert from matrix to SBM, we do
#
#       matrix -->> sparse -->> SparseBlockMatrixR
#
# In theory, this class can be used externally as a useful data structure for storing sparse matrices
#   as an alternative to the Matrix class provided by the Matrix package. Currently, however, the class
#   structure is fairly limited, so there isn't much a reason to do this.
#
#

#' sparse class
#'
#' Low-level representation of sparse matrices
#'
#' An alternative data structure for storing sparse matrices in R using the (row, column, value)
#' format. Internally it is stored as a list with three components, each vectors, that contain
#' the rows / columns / values of the nonzero elements.
#'
#' @param x Various \code{R} objects.
#' @param ... (optional) additional arguments.
#'
#' @docType class
#' @name sparse
NULL

#------------------------------------------------------------------------------#
# is.sparse
#
#' @rdname sparse
#' @export
is.sparse <- function(x){
    inherits(x, "sparse")
} # END IS.SPARSE

#------------------------------------------------------------------------------#
# reIndexC.sparse
#  Re-indexing TO C for sparse objects
#
# #' @describeIn reIndexC C-style re-indexing for \link{sparse} objects.
#' @export
reIndexC.sparse <- function(x){
    if(x$start == 0){
        warning("This object already uses C-style indexing!")
        return(x)
    }

    x$rows <- x$rows - 1
    x$cols <- x$cols - 1
    x$start <- 0

    x
} # END REINDEXC.SPARSE

#------------------------------------------------------------------------------#
# reIndexR.sparse
#  Re-indexing TO R for sparse objects
#
# #' @describeIn reIndexC R-style re-indexing for \link{sparse} objects.
#' @export
reIndexR.sparse <- function(x){
    if(x$start == 1){
        warning("This object already uses R-style indexing!")
        return(x)
    }

    x$rows <- x$rows + 1
    x$cols <- x$cols + 1
    x$start <- 1

    x
} # END REINDEXR.SPARSE

#------------------------------------------------------------------------------#
# sparse.list
#  List constructor
#
#' @export
sparse.list <- function(x, ...){

    if( !is.list(x)){
        stop("Input must be a list!")
    }

    if( length(x) != 5 || names(x) != c("rows", "cols", "vals", "dim", "start") || is.null(names(x))){
        stop("Input is not coercable to an object of type sparse, check list for the following (named) elements: rows, cols, vals, dim, start")
    }

    if( length(unique(lapply(x[1:3], length))) > 1){
        stop("rows / cols / vals elements have different sizes; should all have the same length (pp)!!")
    }

    if(length(x$dim) != 2){
        stop("dim attribute must have length 2!")
    }

    if(x$start != 0 && x$start != 1){
        stop("start attribute must be 0 (C-style) or 1 (R-style)!")
    }

    if(!is.integer(x$rows) || !is.integer(x$cols)){
        stop("rows / cols must both be integers!")
    }

    if(!is.numeric(x$vals)){
        stop("vals must be numeric!")
    }

    structure(x, class = "sparse")
} # END SPARSE.LIST

#------------------------------------------------------------------------------#
# sparse.matrix
#
#' @export
sparse.matrix <- function(x, index = "R", ...){
    if( nrow(x) != ncol(x)) stop("Input matrix must be square!") # 2-7-15: Why does it need to be square?

    if(index != "R" && index != "C") stop("Invalid entry for index parameter: Must be either 'R' or 'C'!")

    pp <- nrow(x)

    nnz <- which(abs(x) > zero_threshold()) - 1
    vals <- double(length(nnz))
    rows <- integer(length(nnz))
    cols <- integer(length(nnz))
    for(k in seq_along(nnz)){
        col <- trunc(nnz[k] / pp)
        row <- nnz[k] - (pp * col)
        vals[k] <- as.vector(x)[nnz[k] + 1]
        rows[k] <- row
        cols[k] <- col
    }

    sp <- sparse.list(list(rows = as.integer(rows), cols = as.integer(cols), vals = as.numeric(vals), dim = c(pp, pp), start = 0))

    if(index == "R"){
        reIndexR(sp)
    } else{
        sp
    }
} # END SPARSE.MATRIX

#------------------------------------------------------------------------------#
# as.sparse.list
#  Convert FROM list TO sparse
#
#' @export
as.sparse.list <- function(x, ...){
    sparse.list(x)
} # END AS.SPARSE.LIST

#------------------------------------------------------------------------------#
# as.sparse.matrix
#  Convert FROM matrix TO sparse
#  By default, return the object using R indexing. If desired, the method can return C-style indexing by setting
#    index = "C".
#' @export
as.sparse.matrix <- function(x, index = "R", ...){
    sparse.matrix(x, index)
} # END AS.SPARSE.MATRIX

#------------------------------------------------------------------------------#
# as.matrix.sparse
#  Convert FROM sparse TO matrix
#
#' @export
as.matrix.sparse <- function(x, ...){

    if( !is.sparse(x)){
        stop("Input must be a sparse object!")
    }

    if(x$start == 0) x <- reIndexR(x) # if indexing starts at 0, adjust to start 1 instead

    m.dim <- x$dim
    m <- matrix(0, nrow = m.dim[1], ncol = m.dim[2])

    for(k in seq_along(x$vals)){
        m[x$rows[k], x$cols[k]] <- x$vals[k]
    }

    attributes(m)$dim <- x$dim
    # attributes(m)$dimnames <- list()
    rownames(m) <- as.character(1:nrow(m))
    colnames(m) <- as.character(1:ncol(m))

    m
} # END AS.MATRIX.SPARSE

#------------------------------------------------------------------------------#
# as.list.sparse
#  Convert FROM sparse TO list
#
#' @export
as.list.sparse <- function(x, ...){
    list(rows = x$rows, cols = x$cols, vals = x$cols, dim = x$dim, start = x$start)
} # END AS.LIST.SPARSE

#------------------------------------------------------------------------------#
# print.sparse
#  Print function for sparse objects
#  By default, format the output as a three-column matrix [cols | rows | vals] ordered by increasing columns.
#    Optionally, set pretty = FALSE to print the sparse object as a list.
#' @export
print.sparse <- function(x, pretty = TRUE, ...){
    if(pretty){
        out <- cbind(x$cols, x$rows, x$vals)
        colnames(out) <- c("cols", "rows", "vals")
        print(out)
    } else{
        print(as.list(x))
    }

} # END PRINT.SPARSE

#------------------------------------------------------------------------------#
# is.zero.sparse
#  Check to see if a sparse object represents the zero matrix
#
#' @export
is.zero.sparse <- function(x){
    check_if_zero <- (length(x$rows) == 0)

    check_if_zero
} # END IS.ZERO.SPARSE

#------------------------------------------------------------------------------#
# .num_edges.sparse
# Internal function for returning the number of edges in a sparse object
#
.num_edges.sparse <- function(x){
    ### Testing only for now
    if(length(which(abs(x$vals) > zero_threshold())) != length(x$rows)){
        stop("Error in .num_edges.sparse! Please check source code.")
    }

    length(which(abs(x$vals) > zero_threshold()))
} # END .NUM_EDGES.SPARSE
