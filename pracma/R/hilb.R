###
### HILB.R  Hilbert matrix
###


hilb <- function(n) {
    if (!is.numeric(n))
        stop("Input argument 'n' must be a numeric scalar.")
    if (n < 0) return(matrix(NA, nrow=0, ncol=0))
    if (length(n) > 1 || ceiling(n) != floor(n)) {
        n <- floor(n[1])
        warning("Size 'n' should be a single integer number.")
    }

    J <- matrix(rep(1:n, each=n), n, n)
    I <- t(J)
    E <- matrix(1, n, n)
    H <- E / (I + J - 1)
    return(H)
}
