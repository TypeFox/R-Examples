###
### $Id: hilb.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Hilbert matrix.
###


##-----------------------------------------------------------------------------
hilb <- function(n) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric", sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1", sQuote("n")))
    }

    n <- floor(n)
    if (n < 0) {
        return(matrix(NA, nrow=0, ncol=0))
    }

    i <- seq_len(n)
    H <- 1 / outer(i-1, i, "+")

    return(H)
}

