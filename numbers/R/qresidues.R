##
##  q r e s i d u e s . R  Quadratic Residues
##


legendre_sym <- function(a, p) {
    stopifnot(length(a) == 1, length(p) == 1)
    if (!isPrime(p))
        stop("Argument 'p' must be a prime number, at least 2.")

    if (p == 2) {
        S <- if (a %% 2 == 0) 0 else 1

    } else if (a %% p == 0) {
        S <- 0

    } else {
        S <- modpower(a, (p-1)/2, p)
        S <- as.numeric(S)
        if (S != 1) S <- -1
    }
    return(S)
}


jacobi_sym <- function(a, n) {
    P <- primeFactors(n)
    S <- 1
    for (p in P) S <- S * legendre_sym(a, p)
    return(S)
}


quadratic_residues <- function(n) {
    stopifnot(length(n) == 1, is.numeric(n))

    if (floor(n) != ceiling(n) || n <= 1)
        stop("Argument 'n' must be a whole number.")
    if (abs(n) > 2^25)
        stop("Argument 'n' may not be larger than 2^25.")

    N <- 0:(abs(n)/2)
    sort(unique(mod(N^2, n)))
}
