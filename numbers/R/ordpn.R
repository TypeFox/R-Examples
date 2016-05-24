ordpn <- function(p, n) {
    # compute the order of p in n!
    stopifnot(is.numeric(p), length(p) == 1,
              is.numeric(n), length(n) == 1)
    if (floor(n) != ceiling(n) || n < 1 ||
        floor(p) != ceiling(p) || p < 1 )
        stop("Arguments 'p' and  'n' must be natural numbers.")
    if (!isPrime(p))
        stop("Argument 'p' must be a prime number.")

    rs = 0
    pk = n/p                # pk = p
    while (pk >= 1) {       # while (n > pk)
        rs = rs + floor(pk) # rs = rs + floor(n/pk)
        pk = pk/p           # pk = p * pk
    }
    return(rs)
}

