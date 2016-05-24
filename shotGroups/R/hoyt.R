#####---------------------------------------------------------------------------
## implement recycling rule for function arguments
#####---------------------------------------------------------------------------

recycle <-
function(...) {
    dots <- list(...)
    maxL <- max(vapply(dots, length, integer(1)))
    lapply(dots, rep, length=maxL)
}

#####---------------------------------------------------------------------------
## Hoyt / Nakagami-q distribution
## correlated bivariate normal distribution rewritten in polar coordinates
## pdf, cdf, and inverse cdf of the distribution of the radius
#####---------------------------------------------------------------------------

## determine parameters for Hoyt distribution
getHoytParam <-
function(x) {
    UseMethod("getHoytParam")
}

## based on data frame with (x,y)-coords
getHoytParam.data.frame <-
function(x) {
    sigma <- cov(getXYmat(x))            # covariance matrix
    x     <- eigen(sigma)$values         # eigenvalues
    NextMethod("getHoytParam")
}

## based on list of covariance matrices
getHoytParam.list <-
function(x) {
    if(!all(vapply(x, is.matrix,  logical(1)))) { stop("x must be a matrix") }
    if(!all(vapply(x, is.numeric, logical(1)))) { stop("x must be numeric") }
    if(!all(vapply(x, dim, integer(2)) == 2L))  { stop("x must be (2 x 2)-matrix") }

    getEV <- function(sigma) {           # eigenvalues from covariance matrix
        if(!isTRUE(all.equal(sigma, t(sigma)))) {
            stop("x must be symmetric")
        }

        lambda <- eigen(sigma)$values
        if(!all(lambda >= -sqrt(.Machine$double.eps) * abs(lambda[1]))) {
            stop("x is numerically not positive definite")
        }
        lambda
    }

    ev    <- lapply(x, getEV)              # eigenvalues for all matrices
    ev1   <- vapply(ev, head, FUN.VALUE=numeric(1), n=1)  # all first eigenvalues
    ev2   <- vapply(ev, tail, FUN.VALUE=numeric(1), n=1)  # all second eigenvalues
    qpar  <- 1/sqrt(((ev1+ev2)/ev2) - 1)   # Hoyt q
    omega <- ev1+ev2                       # Hoyt omega

    return(list(q=qpar, omega=omega))
}

## based on covariance matrix
getHoytParam.matrix <-
function(x) {
    if(any(dim(x) != 2L))           { stop("x must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(x, t(x)))) { stop("x must be symmetric") }

    x <- eigen(x)$values
    NextMethod("getHoytParam")
}

## based on 2-vector of eigenvalues
## not vectorized
getHoytParam.default <-
function(x) {
    if(!is.numeric(x))  { stop("x must be numeric") }
    if(any(x < 0))      { stop("x must be >= 0") }
    if(length(x) != 2L) { stop("x must have length 2") }
    if(!all(x >= -sqrt(.Machine$double.eps) * abs(max(x)))) {
        stop("x is numerically not positive definite")
    }

    x   <- sort(x, decreasing=TRUE)          # largest eigenvalue first
    ev1 <- x[1]
    ev2 <- x[2]

    qpar  <- 1 / sqrt(((ev1+ev2) / ev2) - 1) # Hoyt q
    omega <- ev1+ev2                         # Hoyt omega

    return(list(q=qpar, omega=omega))
}

# determine eigenvalues from Hoyt parameters
getEVfromHoyt <-
function(qpar, omega) {
    nnaQ <- which(!is.na(qpar))
    nnaO <- which(!is.na(omega))
    stopifnot(all(qpar[nnaQ] > 0), all(qpar[nnaQ] < 1), all(omega[nnaO] > 0))

    ev2 <- omega / ((1/qpar^2) + 1)      # 2nd eigenvalue
    ev1 <- omega - ev2                   # 1st eigenvalue

    ## sort each pair of eigenvalues in descending order
    ev1ord <- pmax(ev1, ev2)
    ev2ord <- pmin(ev1, ev2)

    return(list(ev1=ev1ord, ev2=ev2ord))
}

#####---------------------------------------------------------------------------
## pdf Hoyt distribution
## http://reference.wolfram.com/mathematica/ref/HoytDistribution.html
dHoyt <-
function(x, qpar, omega) {
    is.na(x)     <- is.nan(x)                # replace NaN with NA
    is.na(qpar)  <- (qpar  <= 0) | (qpar >= 1) | !is.finite(qpar)
    is.na(omega) <- (omega <= 0) | !is.finite(omega)

    argL  <- recycle(x, qpar, omega)
    x     <- argL[[1]]
    qpar  <- argL[[2]]
    omega <- argL[[3]]

    dens <- numeric(length(x))                 # initialize density to 0
    keep <- which((x >= 0) | !is.finite(x))    # keep non-negative x, NA, -Inf, Inf
    if(length(keep) < 1L) { return(dens) }     # nothing to do

    lfac1 <- log(x[keep]) + log(1 + qpar[keep]^2) - log(qpar[keep]*omega[keep])
    lfac2 <- -x[keep]^2*(1+qpar[keep]^2)^2/(4*qpar[keep]^2*omega[keep])
    bArg  <- (x[keep]^2*(1-qpar[keep]^4)  /(4*qpar[keep]^2*omega[keep]))
    lfac3 <- log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg
    res   <- exp(lfac1+lfac2+lfac3)            # this may be NaN
    dens[keep] <- ifelse(is.nan(res), 0, res)  # if so, set to 0

    return(dens)
}

## equivalent
## Hoyt, RS. 1947. Probability functions for the modulus and angle of the
## normal complex variate. Bell System Technical Journal, 26(2). 318-359.
## Hoyt pdf is for scaled variables with S := 1/sqrt(Su^2+Sv^2), u=U/S, v=V/S
## -> set r to r/S and pdf to pdf/S
# dCNhoyt <- function(r, sigma) {
#     ev <- eigen(sigma)$values
#     b  <- abs(diff(ev)) / sum(ev)
#     S  <- sqrt(sum(ev))
#     r  <- r/S
#
#     fac1 <- (2*r/sqrt(1-b^2)) * exp(-r^2/(1-b^2))
#     bArg <- (b*r^2/(1-b^2))
#     fac2 <- exp(log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg)
#     dens <- fac1*fac2 / S
#
#     return(dens)
# }

## equivalent
## Greenwalt, CR & Shultz, ME. 1968.
## Principles of Error Theory and Cartographic Applications
## ACIC TR-96, Appendix D-3, eq. 3
# dGreenwalt <- function(r, sigma) {
#     ev   <- eigen(sigma)$values
#     fac1 <- 1/prod(sqrt(ev))
#     fac2 <- r*exp(-(r^2/(4*ev[1])) * (1 + (ev[1]/ev[2])))
#     bArg <-        (r^2/(4*ev[1])) * ((ev[1]/ev[2]) - 1)
#     fac3 <- exp(log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg)
#     dens <- fac1*fac2*fac3
#
#     return(dens)
# }

#####---------------------------------------------------------------------------
## generalized Marcum Q-function from non-central chi^2 distribution
## Nuttall, AH. (1975). Some integrals involving the Q-M function.
## IEEE Transactions on Information Theory, 21 (1), 95-96
marcumQ <-
function(a, b, nu, lower.tail=TRUE) {
    pchisq(b^2, df=2*nu, ncp=a^2, lower.tail=lower.tail)
}

#####---------------------------------------------------------------------------
## cdf Hoyt distribution in closed form
## Paris, JF. 2009. Nakagami-q (Hoyt) distribution function with applications.
## Electronics Letters, 45(4). 210-211. Erratum: doi:10.1049/el.2009.0828
pHoyt <-
function(q, qpar, omega, lower.tail=TRUE) {
    is.na(qpar)  <- (qpar  <= 0) | (qpar >= 1) | !is.finite(qpar)
    is.na(omega) <- (omega <= 0) | !is.finite(omega)

    argL  <- recycle(q, qpar, omega)
    q     <- argL[[1]]
    qpar  <- argL[[2]]
    omega <- argL[[3]]

    pp   <- numeric(length(q))               # initialize probabilities to 0
    keep <- which((q >= 0) | !is.finite(q))  # keep non-negative q, NA, NaN, -Inf, Inf

    alphaQ <- (sqrt((1 - qpar[keep]^4))/(2*qpar[keep])) * sqrt((1 + qpar[keep])/(1 - qpar[keep]))
     betaQ <- (sqrt((1 - qpar[keep]^4))/(2*qpar[keep])) * sqrt((1 - qpar[keep])/(1 + qpar[keep]))

    y <- q[keep] / sqrt(omega[keep])
    if(lower.tail) {
        pp[keep] <- marcumQ( betaQ*y, alphaQ*y, nu=1, lower.tail=lower.tail) -
                    marcumQ(alphaQ*y,  betaQ*y, nu=1, lower.tail=lower.tail)

        ## special cases not caught so far
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[keep] <- 1 + marcumQ( betaQ*y, alphaQ*y, nu=1, lower.tail=lower.tail) -
                        marcumQ(alphaQ*y,  betaQ*y, nu=1, lower.tail=lower.tail)

        ## special cases not caught so far
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }

    return(pp)
}

## equivalent
## Hoyt, RS. 1947. Probability functions for the modulus and angle of the
## normal complex variate. Bell System Technical Journal, 26(2). 318-359.
# pCNhoyt <- function(qq, sigma) {
#     ev <- eigen(sigma)$values
#     b  <- abs(diff(ev)) / sum(ev)
#     S  <- sqrt(sum(ev))
#     qq <- qq/S                           # rescale
#
#     intFun <- function(r, b) {
#         fac1 <- r*exp(-(r^2/(1-b^2)))
#         bArg <- (b*r^2/(1-b^2))
#         fac2 <- exp(log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg)
#         res  <- fac1*fac2                # this may be NaN
#         ifelse(is.finite(res), res, 0)   # if so, return 0
#     }
#
#     pp <- (1/sqrt(1-b^2)) * sapply(qq, function(x) 2*integrate(intFun, 0, x, b=b)$value)
#     return(pp)
# }

## equivalent
## Greenwalt, CR & Shultz, ME. 1968.
## Principles of Error Theory and Cartographic Applications
## ACIC TR-96, Appendix D-3, eq3
# pCNgreenwalt <- function(qq, sigma) {
#     intFun <- function(r, ev) {
#         fac1 <- r*exp(-(r^2/(4*ev[1])) * (1 + (ev[1]/ev[2])))
#         ## modified Bessel function of first kind and order 0
#         bArg <-   (r^2/(4*ev[1])) * ((ev[1]/ev[2]) - 1)
#         fac2 <- exp(log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg)
#         res  <- fac1*fac2                        # this may be NaN
#         return(ifelse(is.finite(res), res, 0))   # if so, return 0
#     }
#
#     ev <- eigen(sigma)$values
#     pp <- (1/prod(sqrt(ev))) * sapply(qq, function(x) integrate(intFun, 0, x, ev=ev)$value)
#     return(pp)
# }

## equivalent
## Hoover, WE. 1984. Algorithms For Confidence  Circles, and Ellipses.
## Washington, D.C., National Oceanic and Atmospheric Administration.
## NOAA Technical Report NOS 107 C&GS 3, 1-29. p. 9.
# pCNhoover <- function(qq, sigma) {
#     ev     <- eigen(sigma)$values
#     Hk     <- qq / sqrt(ev[1])
#     Hc     <- sqrt(ev[2] / ev[1])
#     Hbeta  <- 2*Hc / pi
#     Hgamma <- (Hk/(2*Hc))^2
#
#     Hw <- function(phi, Hc) {
#         (Hc^2 - 1)*cos(phi) - (Hc^2 + 1)
#     }
#
#     Hf <- function(phi, Hc, Hgamma) {
#         (exp(Hgamma*Hw(phi, Hc)) - 1) / Hw(phi, Hc)
#     }
#
#     Hbeta * integrate(Hf, 0, pi, Hc=Hc, Hgamma=Hgamma)$value
# }

#####---------------------------------------------------------------------------
## Hoyt quantile function through root finding of cdf
qHoyt <-
function(p, qpar, omega, lower.tail=TRUE, loUp=NULL) {
    is.na(qpar)  <- (qpar <= 0)  | (qpar >= 1) | !is.finite(qpar)
    is.na(omega) <- (omega <= 0) | !is.finite(omega)

    argL  <- recycle(p, qpar, omega)
    p     <- argL[[1]]
    qpar  <- argL[[2]]
    omega <- argL[[3]]

    qq   <- rep(NA_real_, length(p))
    keep <- which((p >= 0) & (p < 1))
    if(length(keep) < 1) { return(qq) }

    if(is.null(loUp)) {                  # no search interval given
        ## use Grubbs chi^2 quantile for setting root finding interval
        ## Grubbs-Liu chi^2 and Hoyt can diverge
        GP <- getGPfromHP(qpar, omega)   # Grubbs parameters
        qGrubbs   <- qChisqGrubbs(p[keep], m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                  varX=GP$varX, l=GP$l, delta=GP$delta,
                                  lower.tail=lower.tail, type="Liu")
        qLo  <- ifelse(p[keep] <= 0.5, 0,         0.25*qGrubbs)
        qUp  <- ifelse(p[keep] <= 0.5, qGrubbs.6, 1.75*qGrubbs)
        loUp <- split(cbind(qLo, qUp), seq_along(p))
    } else {
        if(is.matrix(loUp)) {
            loUp <- split(loUp, seq_len(nrow(loUp)))
        } else if(is.vector(loUp)) {
            loUp <- list(loUp)
        } else if(!is.list(loUp)) {
            stop("loUp must be a list, a matrix, a vector, or missing entirely")
        }
    }

    cdf <- function(x, p, qpar, omega, lower.tail) {
        pHoyt(x, qpar=qpar, omega=omega, lower.tail=lower.tail) - p
    }

    getQ <- function(p, qpar, omega, loUp, lower.tail) {
        tryCatch(uniroot(cdf, interval=loUp, p=p, qpar=qpar, omega=omega,
                         lower.tail=lower.tail)$root,
                 error=function(e) return(NA_real_))
    }

    qq[keep] <- unlist(Map(getQ, p=p[keep], qpar=qpar[keep], omega=omega[keep],
                           loUp=loUp[keep], lower.tail=lower.tail[1]))
    return(qq)
}

#####---------------------------------------------------------------------------
## random numbers from Hoyt distribution
rHoyt <-
function(n, qpar, omega, method=c("eigen", "chol", "cdf"), loUp=NULL) {
    is.na(qpar)  <- (qpar <= 0)  | (qpar >= 1) | !is.finite(qpar)
    is.na(omega) <- (omega <= 0) | !is.finite(omega)
    method <- match.arg(method)

    ## if n is a vector, its length determines number of random variates
    n     <- if(length(n) > 1L) { length(n) } else { n }
    qpar  <- qpar[1]                     # only first shape parameter is used
    omega <- omega[1]                    # only first scale parameter is used

    rn <- if(method == "eigen") {
        lambda <- unlist(getEVfromHoyt(qpar, omega))   # eigenvalues

        ## simulated 2D normal vectors with mean 0
        X  <- matrix(rnorm(n*length(lambda)), nrow=n)  # with identity cov-mat
        xy <- X %*% diag(sqrt(lambda), length(lambda))
        sqrt(rowSums(xy^2))              # distances to center
    } else if(method == "chol") {
        lambda <- getEVfromHoyt(qpar, omega)
        sigma  <- cbind(c(lambda$ev1, 0), c(0, lambda$ev2))
        CF     <- chol(sigma, pivot=TRUE) # Cholesky-factor
        idx    <- order(attr(CF, "pivot"))
        CFord  <- CF[, idx]

        ## simulated 2D normal vectors with mean 0
        xy <- matrix(rnorm(n*ncol(sigma)), nrow=n) %*% CFord
        sqrt(rowSums(xy^2))              # distances to center
    } else if(method == "cdf") {
        ## root finding of pHoyt() given uniform random probabilities:
        ## find x such that F(x) - U = 0
        cdf <- function(x, u, qpar, omega) {
            pHoyt(x, qpar=qpar, omega=omega) - u
        }

        ## find quantile via uniroot() with error handling
        getQ <- function(u, qpar, omega, loUp) {
            tryCatch(uniroot(cdf, interval=loUp, u=u, qpar=qpar, omega=omega)$root,
                     error=function(e) return(NA_real_))
        }

        u <- runif(n)                        # uniform random numbers

        ## determine search interval(s) for uniroot()
        if(is.null(loUp)) {                  # no search interval given
            ## use Grubbs chi^2 quantile for setting root finding interval
            ## Grubbs-Liu chi^2 and Hoyt can diverge
            GP <- getGPfromHP(qpar, omega)   # Grubbs parameters and quantiles
            qGrubbs   <- qChisqGrubbs(u, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qGrubbs.6 <- qChisqGrubbs(0.6, m=GP$m, v=GP$v, muX=GP$muX,
                                      varX=GP$varX, l=GP$l, delta=GP$delta, type="Liu")
            qLo  <- ifelse(u <= 0.5, 0,         0.25*qGrubbs)
            qUp  <- ifelse(u <= 0.5, qGrubbs.6, 1.75*qGrubbs)
            loUp <- split(cbind(qLo, qUp), seq_along(u))
        } else {
            if(is.matrix(loUp)) {
                loUp <- split(loUp, seq_len(nrow(loUp)))
            } else if(is.vector(loUp)) {
                loUp <- list(loUp)
            } else if(!is.list(loUp)) {
                stop("loUp must be a list, a matrix, a vector, or missing entirely")
            }
        }

        unlist(Map(getQ, u=u, qpar=qpar, omega=omega, loUp=loUp))
    }

    return(rn)
}
