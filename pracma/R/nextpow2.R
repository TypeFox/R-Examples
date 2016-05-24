###
### NEXTPOW2.R  Next higher power of 2
###

nextpow2 <- function(x) {
    if (is.null(x) || length(x) == 0) return(c())
    if (!is.numeric(x) && !is.complex(x))
        stop("Argument 'x' must be a numeric/complex vector/matrix.")

    x[x == 0] <- 1
    return(ceiling(log2(abs(x))))
}
