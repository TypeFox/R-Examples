###
### $Id: factors.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Factorize natural number.
###


##-----------------------------------------------------------------------------
factors <- function(n) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("n")))
    } else if (!(n >= 0)) {
        stop(sprintf("argument %s must be a nonnegative quantity",
                     sQuote("n")))
    }

    n <- floor(n)
    if (n < 4) {
        return(n)
    }

    p <- matlab::primes(sqrt(n))
    # d <- matlab::find(matlab::rem(n, p) == 0)
    d <- which(n %% p == 0)
    if (length(d) == 0) {
        return(n)  # n is prime
    }

    f <- c()
    for (q in p[d]) {
        while (n %% q == 0) {
            f <- c(f, q)
            n <- n/q
        }
    }
    if (n > 1) {
        f <- c(f, n)
    }

    return(f)
}

