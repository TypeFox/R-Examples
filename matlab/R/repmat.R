###
### $Id: repmat.R 55 2014-02-06 16:41:28Z plroebuck $
###
### Replicate and tile an array.
###


##-----------------------------------------------------------------------------
repmat <- function(A, ...) {
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

    B <- switch(EXPR = mode(A),
                logical =,
                complex =,
                numeric = {
                    if (all(dims == 1)) {
                        ## nothing to do in this case
                        A
                    } else if (dims[length(dims)] == 1) {
                        ## avoid dropping single dimension before replication
                        t(kronecker(array(1, rev(dims)), A))
                    } else {
                        kronecker(array(1, dims), A)
                    }
                },
                character = {
                    ## convert unique strings to factors
                    fA <- factor(A, levels = unique(A))
                    ## recursive call using integer version of A
                    iA.mat <- Recall(unclass(fA), dims)
                    ## replace integers with strings according to factor values
                    saved.dim <- dim(iA.mat)
                    cA.mat <- sapply(seq(along = iA.mat),
                                     function(i, A, fac) {
                                         A[i] <- levels(fac)[A[i]]
                                     },
                                     iA.mat, fA)
                    dim(cA.mat) <- saved.dim
                    cA.mat
                },
                NULL)

    if (is.null(B)) {
        stop(sprintf("argument %s must be one of [%s]",
                     sQuote("A"),
                     paste(c("numeric",
                             "logical",
                             "complex",
                             "character"), collapse  = "|")))
    }

    return(B)
}

