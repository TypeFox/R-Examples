###
### f a c t o r s . R  Factorize natural number
###


factors <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || n != round(n) || n < 1)
        stop("Argument 'n' must be a nonnegative integer.")
    if (n >= 2^53)
        stop("Argument 'n' must not be larger than 2^53-1.")
    if (n < 4) return(n)

    f <- c()
    p <- primes(sqrt(n))
    d <- which(n %% p == 0)
    if (length(d) == 0) return(n)  # n is prime

    for (q in p[d]) {
        while (n %% q == 0) {
            f <- c(f, q)
            n <- n/q
        }
    }
    if (n > 1) f <- c(f, n)
    return(f)
}
