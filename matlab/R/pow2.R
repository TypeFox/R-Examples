###
### $Id: pow2.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Raise 2 to some power.
###


##-----------------------------------------------------------------------------
pow2 <- function(f, e) {
    if (!(is.numeric(f) || is.complex(f))) {
        stop(sprintf("argument %s must be numeric or complex",
                     sQuote("f")))
    }

    if (missing(e)) {
        e <- f
        f <- rep(1, length(e))
    } else {
        if (!(is.numeric(e) || is.complex(e))) {
            stop(sprintf("Argument %s must be numeric or complex",
                         sQuote("e")))
        }

        if (is.complex(f) || is.complex(e)) {
            f <- Re(f)
            e <- Re(e)
            warning("imaginary part of complex arguments ignored")
        }
    }

    if (length(f) != length(e)) {
        stop(sprintf("Arguments %s and %s must be of same length",
                     sQuote("f"),
                     sQuote("e")))
    }

    return(f * 2^e)
}

