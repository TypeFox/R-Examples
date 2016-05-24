###
###  p r i m e s . R  Prime numbers
###


##  Sieve of Erathostenes
primeSieve <- function(n) {
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 1)
        stop("Argument 'n' must be an integer number greater or equal 1.")
    if (n > 2^53 - 1)
        stop("Argument 'n' must be smaller than 2^53 - 1.")

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


Primes <- function(n1 = 1, n2 = NULL) {
    if (is.null(n2))
        return(primeSieve(n1))

    if (!is.numeric(n1) || length(n1) != 1 || floor(n1) != ceiling(n1) || n1 <= 0 ||
        !is.numeric(n2) || length(n2) != 1 || floor(n2) != ceiling(n2) || n2 <= 0 )
        stop("Arguments 'n1' and 'n2' must be integers.")

    if (n2 > 2^53 - 1)  stop("Upper bound 'n2' must be smaller than 2^53-1.")
    if (n1 > n2)        stop("Upper bound must be greater than lower bound.")
   
    if (n2 <= 1000) {
        P <- primeSieve(n2)
        return(P[P >= n1])
    }

    myPrimes <- primeSieve(floor(sqrt(n2)))

    N <- seq.int(n1, n2)
    n <- length(N)
    A <- numeric(n)
    if (n1 == 1) A[1] <- -1

    for (p in myPrimes) {
        r <- n1 %% p                                    # rest modulo p
        if (r == 0) { i <- 1 } else { i <- p - r + 1 }  # find next divisible by p
        if (i <= n && N[i] == p) { i <- i + p }         # if it is p itself, skip
        while (i <= n) { A[i] <- 1; i <- i + p }        # mark those divisible by p
    }
    return(N[A == 0])
}


twinPrimes <- function(n1, n2) {
    P <- Primes(n1, n2)
    twins <- which(diff(P) == 2)
    cbind(P[twins], P[twins+1])
}


nextPrime <-function(n) {
    if (n <= 1)  n <- 1  else  n <- floor(n)
    n <- n + 1

    # m <- 2*n  # Bertrands law
    d1 <- max(3, round(log(n)))
    P  <- Primes(n, n + d1)

    while(length(P) == 0) {
        n <- n + d1 + 1
        P  <- Primes(n, n + d1)
    }
    return( as.numeric(min(P)) )
}


previousPrime <-function(n) {
    if (n <= 2) return(c())
    if (floor(n) == ceiling(n))  n <- n - 1
    else n <- floor(n)

    if (n <= 10) {
        P <- c(2, 3, 5, 7)
        return(max(P[P <= n]))
    }

    # m <- 2*n  # Bertrands law
    d1 <- max(3, round(log(n)))
    P  <- Primes(n - d1, n)

    while(length(P) == 0 || n - d1 < 3) {
        n <- n - d1 - 1
        P  <- Primes(n - d1, n)
    }
    return( as.numeric(max(P)) )
}


atkin_sieve <- function(n) {
    stopifnot(length(n) == 1, floor(n) == ceiling(n), n >= 1)

    if (n <= 4)
        return(which(c(FALSE, TRUE, TRUE, FALSE)[1:n]))

    sieve <- vector(mode = "logical", length = n)
    sqrtn <- floor(sqrt(n))
    for (x in 1:sqrtn) {
        for (y in 1:sqrtn) {
            m <- 4*x^2 + y^2
            if (m <= n && (m %% 12 == 1 || m %% 12 == 5))
                sieve[m] <- !sieve[m]
  
            m <- 3*x^2 + y^2
            if (m <= n && m %% 12 == 7)
                sieve[m] <- !sieve[m]
  
            m <- 3*x^2 - y^2
            if (x > y && m <= n && m %% 12 == 11)
                sieve[m] <- !sieve[m]
        }
    }
    if (n >= 25) {
        for (x in 5:sqrtn) {
            if (sieve[x])
                sieve[seq.int(x^2, n, by = x^2)] <- FALSE
        }
    }

    sieve[c(2, 3)] <- TRUE
    return(which(sieve))
}
