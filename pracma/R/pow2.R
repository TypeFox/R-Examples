###
### POW2.R  Raise 2 to some power
###

pow2 <- function(f, e) {
    if (!is.numeric(f) && !is.complex(f))
        stop("Argument 'f' must be numeric or complex.")
    if (missing(e)) {
        e <- f
        f <- rep(1, length(e))
    } else {
        if (!is.numeric(e) && !is.complex(e))
            stop("Argument 'e' must be numeric or complex.")
        if (is.complex(f) || is.complex(e)) {
            f <- Re(f)
            e <- Re(e)
            warning("Imaginary part of arguments ignored.")
        }
    }
    if (length(f) != length(e))
        stop("Arguments 'e' and 'f' must be of same length.")

    return(f * 2^e)
}
