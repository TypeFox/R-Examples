###
### f a c t o r i z e . R  Factorize natural number
###


primeFactors <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || n != round(n) || n < 1) {
        warning("Argument 'n' must be a nonnegative integer.")
        return(NULL)
    }
    if (n < 4) return(n) 

    if (n <= 2^53 - 1) {
        f <- c()
        p <- Primes(floor(sqrt(n)))
        d <- which(n %% p == 0)
        if (length(d) == 0) return(n)  # n is prime
    
        for (q in p[d]) {
            while (n %% q == 0) {
                f <- c(f, q)
                n <- n/q
            }
        }
        if (n > 1) f <- c(f, n)

    } else {
        warning("Argument 'n' too big: use 'gmp::factorize()' for this.")
        f <- NA
    }

    return(f)
}


radical <- function(n) {
    prod(unique(primeFactors(n)))
}

