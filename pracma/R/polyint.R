###
### POLYINT.R  Polynom
###

polyint <- function(p, k=0) {
    if (length(p) == 0) return(c())
    if (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex"))
        stop("Argument 'p' must be a real or complex vector.")
    if (!is.vector(k, mode="numeric") && !is.vector(k, mode="complex"))
        stop("Argument 'k' must be a real or complex vector")

    return( c(p / (length(p):1), k) )
}
