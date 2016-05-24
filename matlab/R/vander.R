###
### $Id: vander.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Returns the Vandermonde matrix.
###


##-----------------------------------------------------------------------------
vander <- function(v) {
    if (!is.vector(v)) {
        stop(sprintf("argument %s must be vector", sQuote("v")))
    } else if (!(is.numeric(v) || is.complex(v))) {
        stop(sprintf("argument %s must be numeric or complex", sQuote("v")))
    }

    n <- length(v)
    A <- if (n == 0) {
             matrix(as.numeric(NA), nrow=0, ncol=0)
         } else {
             outer(v, seq(n-1, 0), "^")
         }
    return(A)
}

