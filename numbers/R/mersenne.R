##
##  m e r s e n n e . R  Mersenne Numbers
##


mersenne <- function(p) {
    stopifnot(is.numeric(p), length(p) == 1)
    if (!isNatural(p) || !isPrime(p))
        stop("Argument 'p' must be a prime number for 2^p-1 to be prime.")
    if (p == 2)
        return(TRUE)

    if (!requireNamespace("gmp", quietly = TRUE)) {
        stop("Package 'gmp' needed: Please install separately.", call. = FALSE)
    }


    z2 <- gmp::as.bigz(2)
    z4 <- z2 * z2
    zp <- gmp::as.bigz(p)
    zm <- z2^zp - 1                     # candidate Mersenne prime
    S  <- rep(z4, p - 1)

    for (n in 1:(p-2))
        S[n+1] <- gmp::mod.bigz(S[n]*S[n] - z2, zm)

    if (S[p-1] == 0) tf <- TRUE         # Lucas-Lehmer test
    else             tf <- FALSE
    return(tf)
}
