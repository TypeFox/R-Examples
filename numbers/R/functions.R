##
##  n u m b e r s . R  Functions from Number Theory
##


eulersPhi <- function(n) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")
    if (n == 1) return(1)

    m <- n
    for (p in unique(primeFactors(n)))
        m <- m * (1 - 1/p)
    return(round(m))
}


moebius <- function(n) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")

    R <- rle(primeFactors(n))
    if (n == 1) {
        r <- 1
    } else if (max(R$lengths) > 1) {
        r <- 0
    } else {
        r <- (-1)^length(R$values)
    }
    return(r)
}


mertens <- function(n) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")

    sum(sapply(1:n, moebius))
}


Sigma <- function(n, k = 1, proper = FALSE) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")
    if (!is.numeric(k) || length(k) != 1)
        stop("Argument 'k' must be a numeric scalar.")

    if (n == 1) return(if (proper) 0 else 1)
    R <- rle(primeFactors(n))
    P <- 1
    for (i in 1:length(R$values)) {
        ri <- R$values[i]
        ai <- R$lengths[i]
        P <- P * sum(rep(ri, ai+1)^seq(0, ai*k, length=ai+1))
    }
    if (proper) P <- P - n^k
    return(P)
}


tau <- function(n) {  # Ramanujan's tau function
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")

    if (n == 0) s <- 0
    else if (n == 1) s <- 1
    else {
        s <- 0
        for (k in 1:(n-1)) {
            s <- s + Sigma(k, 5)*Sigma(n-k, 5)
        }
        s <- 65/756 * Sigma(n, 11) + 691/756 * Sigma(n, 5) - 691/3 * s
    }
    return(s)
}


omega <- function(n) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")

    if (n == 1) 0
    else length(unique(primeFactors(n)))
}


Omega <- function(n) {
    if (!isNatural(n))
        stop("Argument 'n' must be a single positive integers.")

    if (n == 1) 0
    else sum(rle(primeFactors(n))$length)
}
