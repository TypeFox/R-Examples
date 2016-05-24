#' Create a cell array (a matrix of lists) 
#' 
#' Create a cell array as in MATLAB.
#' 
#' @param ... Dimensions for the cell array.
#' 
#' @return A list of empty lists.
#' 
#' @note This is a modified version of the same function in \code{\link[matlab]{cell}}
#' @family MATLAB
#' 
#' @example tests/testthat/examples_fcn_doc/examples_cell.R
#' 
#' @export
## Function written to match MATLAB function
## Author: Andrew Hooker

cell <- function (...) 
{
    nargs <- length(dots <- list(...))
    dims <- as.integer(if (nargs == 1) {
      ##dims <- as.integer(if (nargs == 1 && matlab:::is.size_t(dots[[1]])) {
      dots[[1]]
    } else {
        unlist(dots)
    })
    if (length(dims) == 1) {
        dims[2] <- dims[1]
    }
    if (!(length(dims) > 1)) {
        stop("dimensions must be of length greater than 1")
    }
    else if (!(all(dims > 0))) {
        ##stop("dimensions must be a positive quantity")
        L <- list()
        dim(L) <- dims
        return(L)
    }
    createCellArray <- function(dims) {
        L <- vector("list", prod(dims))
        dim(L) <- dims
        for (i in seq(along = L)) {
            L[[i]] <- numeric(0)
        }
        return(L)
    }
    return(createCellArray(dims))
}
