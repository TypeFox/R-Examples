dTriangular <- function(a, b, c, d, e){

    ## computes the density at all x in {a, ..., d}
    x <- seq(a, d, by = 1)
    n <- length(x)
    f <- rep(NA, n)

    # left part
    ind1 <- (1:n)[(x >= a) & (x < b)]
    f[ind1] <- c / (b - a) * (x[ind1]  - a)

    # right part
    ind2 <- (1:n)[(x >= b) & (x <= d)]
    f[ind2] <- (e - c) / (d - b) * (x[ind2] - b) + c

    ## f is the log-density
    dens <- exp(f)
    dens <- dens / sum(dens)
    return(dens)
}
