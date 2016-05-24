###
### $Id: primes.R 48 2014-02-05 20:50:54Z plroebuck $
###
### Generate list of prime numbers.
###


##-----------------------------------------------------------------------------
primes <- function(n) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("n")))
    }

    n <- floor(n)
    if (n < 2) {
        return(c())
    }

    p <- seq(1, n, by=2)
    p[1] <- 2
    q <- length(p)

    if (n >= 9) {
        for (k in seq(3, sqrt(n), by=2)) {
            if (p[(k+1)/2] != 0) {
                p[seq((k*k+1)/2, q, by=k)] <- 0
            }    
        }    
    }

    p[p > 0]
}


##-----------------------------------------------------------------------------
## Alternative version...
eratosthenes_sieve <- function(n) {
    if (!is.numeric(n)) {
        stop(sprintf("argument %s must be numeric",
                     sQuote("n")))
    } else if (!(length(n) == 1)) {
        stop(sprintf("argument %s must be of length 1",
                     sQuote("n")))
    }

    n <- floor(n)
    if (n < 2) {
        return(c())    ## :TBD: integer(0) instead?
    }

    ## :HACK: have issues with seq(by=) output for [2 <= n <= 5]
    switch(EXPR=as.character(n),
           "2"=return(as.integer(2)),
           "3"=,
           "4"=return(as.integer(c(2, 3))),
           "5"=return(as.integer(c(2, 3, 5))))

    ## Create a candidate list within which non-primes will be
    ## marked as 0; only candidates below sqrt(n) need be checked. 
    candidates <- seq(from=1, to=n)    
    fin <- floor(sqrt(n))

    ## Loop over the candidates, marking out each multiple
    for (i in seq(from=2, to=fin+1)) {
        if (candidates[i]) {
            candidates[seq(from=2*i, to=n, by=i)] <- 0
        }
    }
    
    ## Filter out non-primes
    candidates[candidates > 1]
}

