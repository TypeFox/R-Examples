## Maxwell-Boltzmann distribution
## http://de.scribd.com/doc/23073705/Maxwell-Distribution-Rev3

#####---------------------------------------------------------------------------
## Maxwell-Boltzmann distribution
dMaxwell <-
function(x, sigma=1) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(x, sigma)
    x     <- args[[1]]
    sigma <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1L) { return(dens) }

    dens[keep] <- sqrt(2/pi) * (x^2/sigma[keep]^3)*exp(-x[keep]^2/(2*sigma[keep]^2))
    return(dens)
}

## Maxwell-Boltzmann cdf
pMaxwell <-
function(q, sigma=1, lower.tail=TRUE) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(q, sigma)
    q     <- args[[1]]
    sigma <- args[[2]]

    pp   <- numeric(length(q))
    keep <- which((q >= 0) | !is.finite(q))

    if(lower.tail) {
        pp[keep] <- 2*pnorm(q[keep]/sigma[keep]) - 1 -
            sqrt(2/pi) * (q[keep]/sigma[keep])*exp(-q[keep]^2/(2*sigma[keep]^2))

        ## some special values not caught before
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[keep] <- 1 - (2*pnorm(q[keep]/sigma[keep]) - 1 -
            sqrt(2/pi) * (q[keep]/sigma[keep])*exp(-q[keep]^2/(2*sigma[keep]^2)))

        ## some special values not caught before
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }

    return(pp)
}

## Maxwell-Boltzmann quantile function
qMaxwell <-
function(p, sigma=1, lower.tail=TRUE) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(p, sigma)
    p     <- args[[1]]
    sigma <- args[[2]]

    keep <- which((p >= 0) & (p < 1))
    qq   <- rep(NA_real_, length(p))
    if(length(keep) < 1L) { return(qq) }

    qq[keep] <- sqrt(qgamma(p[keep], shape=3/2, scale=2*sigma[keep]^2,
                            lower.tail=lower.tail))

    return(qq)
}

## random deviates
rMaxwell <-
function(n, sigma=1) {
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    ## simulate coords separately instead of matrix(rnorm(3*n), ncol=3)
    ## for correct cyclic replication of sigma
    xyz <- replicate(3, rnorm(n, mean=0, sd=sigma))
    sqrt(rowSums(xyz^2))
}
