rdirichlet <- function (n, shape) {
    if (length(shape)<2) stop("`length(shape)' must be an integer greater than 1")
    if (length(n)>1) n <- length(n)
    if (length(shape)==2) { # reduced to generate beta variates
        warning("'shape' is of length 2: reduced to beta variates")
        return(rbeta(n,shape[1],shape[2]))
    } else {
        out <- .Call("rdirichlet", n, shape)
        return(out)
    }
}
