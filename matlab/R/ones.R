###
### $Id: ones.R 55 2014-02-06 16:41:28Z plroebuck $
###
### Create a matrix of all ones.
###


##-----------------------------------------------------------------------------
ones <- function(...) {
    nargs <- length(dots <- list(...))
    dims <- as.integer(if (nargs == 1 && is.size_t(dots[[1]])) {
                           dots[[1]]
                       } else {
                           unlist(dots)
                       })

    if (length(dims) == 1) {
        dims[2] <- dims[1]
    }

    if (!(length(dims) > 1)) {
        stop("dimensions must be of length greater than 1")
    } else if (!(all(dims > 0))) {
        stop("dimensions must be a positive quantity")
    }

    return(array(1, dims))
}

