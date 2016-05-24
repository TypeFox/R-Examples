###
### POLYVAL.R  Polynom
###

polyval <- function(p, x) {
    if (length(x) == 0) return(c())
    if (length(p) == 0) return(0 * x)
    if (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex"))
        stop("Argument 'p' must be a real or complex vector.")
    if (!is.vector(x) && !is.matrix(x))
        stop("Argument 'x' must be a real or complex matrix.")

    n <- length(p)
    y <- outer(x[1:length(x)], (n-1):0, "^") %*% p
    dim(y) <- dim(x)
    return(y)
}
