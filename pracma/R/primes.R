###
###  p r i m e s . R  Prime numbers
###


primes <- function(n) {
    if (!is.numeric(n) || length(n) != 1)
        stop("Argument 'n' must be a numeric scalar.")
    n <- floor(n)
    if (n < 2) return(c())
    p <- seq(1, n, by=2)
    q <- length(p)
    p[1] <- 2
    if (n >= 9) {
        for (k in seq(3, sqrt(n), by=2)) {
            if (p[(k+1)/2] != 0)
                p[seq((k*k+1)/2, q, by=k)] <- 0
        }    
    }
    p[p > 0]
}
