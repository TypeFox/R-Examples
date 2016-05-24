###
### $Id: nextpow2.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Next higher power of 2.
###


##-----------------------------------------------------------------------------
nextpow2 <- function(x) {
    if (!(is.numeric(x) || is.complex(x))) {
        stop(sprintf("argument %s must be numeric or complex",
                     sQuote('x')))
    }

    if (length(x) == 0) {
        return(numeric(0))
    }

    x[x == 0] <- 1
    return(ceiling(log2(abs(x))))
}

