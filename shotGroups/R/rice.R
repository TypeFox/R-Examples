#####---------------------------------------------------------------------------
## Rice distribution
## http://reference.wolfram.com/language/ref/RiceDistribution.html
## uncorrelated bivariate normal distribution with equal variances
## rewritten in polar coordinates
## pdf, cdf, and inverse cdf of the distribution of the radius around an
## offset center: nu is the offset (distance to origin), sigma the scale
#####---------------------------------------------------------------------------

## estimate Rice parameters nu, sigma, MR, RSD, from set of 2D coordinates
getRiceParam <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    UseMethod("getRiceParam")
}

getRiceParam.data.frame <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    xy <- getXYmat(xy)
    NextMethod("getRiceParam")
}

getRiceParam.default <-
function(xy, level=0.95, doRob=FALSE, type=c("LiZhangDai", "MOM")) {
    xy <- as.matrix(xy)
    p  <- ncol(xy)
    if(!is.numeric(xy)) { stop("xy must be numeric") }
    if(!(p %in% (1:3))) { stop("x must be (n x 1/2/3)-matrix") }

    type <- match.arg(type)

    ## check if we can do robust estimation if so required
    haveRob <- if(nrow(xy) < 4L) {
        if(doRob) {
            warning("We need >= 4 points for robust estimations")
        }
        FALSE
    } else {
        rob <- robustbase::covMcd(xy, cor=FALSE)
        TRUE
    }                                    # if(nrow(xy) < 4L)

    ctr <- if(doRob && haveRob) {        # estimated center
        rob$center                       # robust estimate
    } else {
        colMeans(xy)                     # coord-wise mean
    }                                    # if(doRob && haveRob)

    ## get sigma estimate from Rayleigh distribution
    sigmaHat <- getRayParam(xy=xy, level=level, doRob=doRob)$sigma

    ## estimate nu^2 -> but E(xBar'xBar) = mu'mu + (p/N)*sigma^2
    N    <- nrow(xy)
    bias <- (p/N)*sigmaHat["sigma"]^2

    ## conventional estimator: xBar'xBar - (p/N)*sigmaHat^2
    ce <- sum(ctr^2) - bias

    nuSqHat <- if(type == "MOM") {
        ## set xBar'xBar - (p/N)*sigmaHat^2 to 0 when (p/N)*sigmaHat^2 > xBar'xBar
        max(ce, 0)
    } else if(type == "LiZhangDai") {
        ## ML -> bad for low SNR -> instead: Li, Zhang & Dai, 2009
        max(ce, (1/(bias+1)) * sum(ctr^2))
    }

    ## c4 correction for negative bias due to taking square root (concave)
    nuHat <- (1/c4(p*N+1))*sqrt(nuSqHat)

    ## radial mean and sd
    MSD <- getMSDfromRice(nu=nuHat, sigma=sigmaHat["sigma"])

    return(list(nu=setNames(nuHat,    NULL),
             sigma=sigmaHat,
                MR=setNames(MSD$mean, NULL),
               RSD=setNames(MSD$sd,   NULL)))
}

## Laguerre half polynomial
LaguerreHalf <-
function(x) {
    a   <- -x/2
    bI0 <- exp(log(besselI(a, nu=0, expon.scaled=TRUE)) + a)
    bI1 <- exp(log(besselI(a, nu=1, expon.scaled=TRUE)) + a)
    exp(x/2) * ((1-x)*bI0 - x*bI1)
}

doubleFactorial <-
function(x, log=FALSE) {
    y   <- (x + 1)/2
    lDF <- lgamma(2*y) - (lgamma(y) + (y-1) * log(2))
    if(log) {
        lDF
    } else {
        exp(lDF)
    }
}

## mean radial error and sd of radial error from Rice parameters nu, sigma
getMSDfromRice <-
function(nu, sigma) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(nu, sigma)
    nu    <- argL[[1]]
    sigma <- argL[[2]]

    L05   <- LaguerreHalf(-0.5 * nu^2 / sigma^2)
    rMean <- sigma * sqrt(pi/2) * L05
    rVar  <- 2*sigma^2 + nu^2 - (pi * sigma^2 / 2) * L05^2

    ## for large signal-to-noise ratios, use approximation (Foi, 2011)
    s2nr <- nu/sigma                           # signal to noise ratio
    rMean[which(s2nr > 52)] <- nu + sigma^2/(2*nu)
    rVar[ which(s2nr > 52)] <- sigma^2 - sigma^4/(2*nu^2)

    return(list(mean=rMean, sd=sqrt(rVar)))
}

#####---------------------------------------------------------------------------
## pdf of Rice distribution
dRice <-
function(x, nu, sigma) {
    is.na(x)     <- is.nan(x)                  # replace NaN with NA
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(x, nu, sigma)
    x     <- argL[[1]]
    nu    <- argL[[2]]
    sigma <- argL[[3]]

    dens <- numeric(length(x))                 # initialize density to 0
    keep <- which((x >= 0) | !is.finite(x))    # keep non-negative x, NA, -Inf, Inf
    if(length(keep) < 1L) { return(dens) }     # nothing to do

    lfac1 <- log(x[keep]) - 2*log(sigma[keep])
    lfac2 <-   -(x[keep]^2 + nu[keep]^2) / (2*sigma[keep]^2)
    bArg  <- abs(x[keep]   * nu[keep]   /     sigma[keep]^2)
    lfac3 <- log(besselI(bArg, nu=0, expon.scaled=TRUE)) + bArg
    res   <- exp(lfac1+lfac2+lfac3)            # this may be NaN
    dens[keep] <- ifelse(is.nan(res), 0, res)  # if so, set to 0

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                           # signal to noise ratio
    rMSD <- getMSDfromRice(nu, sigma)          # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    dens[keepS2NR24] <- dnorm(x[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep])
    dens[keepS2NR52] <- dnorm(x[keep], mean=nu[keep],        sd=sigma[keep])

    return(dens)
}

#####---------------------------------------------------------------------------
## cdf Rice distribution
pRice <-
function(q, nu, sigma, lower.tail=TRUE) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    argL  <- recycle(q, nu, sigma)
    q     <- argL[[1]]
    nu    <- argL[[2]]
    sigma <- argL[[3]]

    pp   <- numeric(length(q))               # initialize probabilities to 0
    keep <- which((q >= 0) | !is.finite(q))  # keep non-negative q, NA, NaN, -Inf, Inf

    aQ <- nu/sigma
    bQ <-  q/sigma
    pp[keep] <- marcumQ(aQ[keep], bQ[keep], nu=1, lower.tail=lower.tail)

    ## special cases not caught so far
    if(lower.tail) {
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                         # signal to noise ratio
    rMSD <- getMSDfromRice(nu, sigma)        # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    pp[keepS2NR24] <- pnorm(q[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep], lower.tail=lower.tail)
    pp[keepS2NR52] <- pnorm(q[keep], mean=nu[keep],        sd=sigma[keep],   lower.tail=lower.tail)

    return(pp)
}

#####---------------------------------------------------------------------------
## Rice quantile function
qRice <-
function(p, nu, sigma, lower.tail=TRUE) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)

    args  <- recycle(p, nu, sigma)
    p     <- args[[1]]
    nu    <- args[[2]]
    sigma <- args[[3]]

    keep <- which((p >= 0) & (p < 1))
    qq   <- sqrt(qchisq(p, df=2, ncp=(nu/sigma)^2)) * sigma # Marcum Q-function

    ## for large signal-to-noise ratios, use normal approximation
    s2nr <- nu/sigma                         # signal to noise ratio
    rMSD <- getMSDfromRice(nu, sigma)        # M, SD of radial error

    keepS2NR24 <- keep[keep %in% which(s2nr > 24)]
    keepS2NR52 <- keep[keep %in% which(s2nr > 52)]
    qq[keepS2NR24] <- qnorm(p[keep], mean=rMSD$mean[keep], sd=rMSD$sd[keep], lower.tail=lower.tail)
    qq[keepS2NR52] <- qnorm(p[keep], mean=nu[keep],        sd=sigma[keep],   lower.tail=lower.tail)

    return(qq)
}

#####---------------------------------------------------------------------------
## random numbers from Rice distribution
rRice <-
function(n, nu, sigma, method=c("eigen", "chol", "cdf")) {
    is.na(nu)    <- (nu < 0)     | !is.finite(nu)
    is.na(sigma) <- (sigma <= 0) | !is.finite(sigma)
    method <- match.arg(method)

    ## if n is a vector, its length determines number of random variates
    n     <- if(length(n) > 1) { length(n) } else { n }
    nu    <- nu[1]                      # only first shape parameter is used
    sigma <- sigma[1]                   # only first scale parameter is used

    rn <- if(method == "eigen") {
        ## simulated 2D normal vectors with mean 0
        X  <- matrix(rnorm(n*2), nrow=n)    # with identity cov-mat
        xy <- X %*% diag(rep(sigma, 2))     # with cov-mat according to sigma

        ## move to mean nu
        xyMove <- sweep(xy, 2, nu, FUN="+")
        sqrt(rowSums(xyMove^2))             # distances to center
    } else if(method == "chol") {
        covMat <- diag(rep(sigma^2, 2))
        CF     <- chol(covMat, pivot=TRUE)  # Cholesky-factor
        idx    <- order(attr(CF, "pivot"))
        CFord  <- CF[, idx]

        ## simulated 2D normal vectors with mean 0 and cov-mat according to sigma
        xy <- matrix(rnorm(n*2), nrow=n) %*% CFord

        ## move to mean nu
        xyMove <- sweep(xy, 2, nu, FUN="+")
        sqrt(rowSums(xyMove^2))          # distances to center
    } else if(method == "cdf") {
        u <- runif(n)                    # uniform random numbers
        sqrt(qchisq(u, df=2, ncp=(nu/sigma)^2)) * sigma
    }

    return(rn)
}
