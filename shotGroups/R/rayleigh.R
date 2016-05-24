#####---------------------------------------------------------------------------
## c4(n) correction factor for taking the square root of the variance of a
## normally distributed random variable with n observations -> E(s) = c4*sigma
## http://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
## http://ballistipedia.com/index.php?title=Measuring_Precision#Estimating_.CF.83
## always <= 1
c4 <-
function(n) {
    stopifnot(all(n >= 2))

    ## use exp(lgamma()) because gamma() will be infinite for large N
    fac <- sqrt(2/(n-1)) * exp(lgamma(n/2) - lgamma((n-1)/2))
    ifelse((fac > 1) | is.infinite(fac), 1, fac)
}

## estimate of Rayleigh sigma parameter
getRayParam <-
function(xy, level=0.95, mu, doRob=FALSE) {
    UseMethod("getRayParam")
}

getRayParam.data.frame <-
function(xy, level=0.95, mu, doRob=FALSE) {
    xy <- getXYmat(xy)
    NextMethod("getRayParam")
}

getRayParam.default <-
function(xy, level=0.95, mu, doRob=FALSE) {
    xy <- as.matrix(xy)
    if(!is.numeric(xy))    { stop("xy must be numeric") }
    if(!is.numeric(level)) { stop("level must be numeric") }
    if(level <= 0)         { stop("level must be > 0") }

    ## check if CI level is given in percent
    if(level >= 1) {
        while(level >= 1) { level <- level / 100 }
        warning(c("level must be in (0,1) and was set to ", level))
    }

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

    N     <- nrow(xy)
    p     <- ncol(xy)
    alpha <- 1-level

    if(!missing(mu)) {                   # Singh C1 - true mean is known
        xyCtr  <- sweep(xy, 2, mu, "-")  # center with true mean
        rSqSum <- if(doRob && haveRob) { # sum squared radii centered data
            ## = N*trace of uncorrected covariance matrix -> N*trace((N-1)/N)*cov(xy)
            (N-1)*sum(diag(rob$cov))
        } else {
            sum(xyCtr^2)                 # sum squared radii
        }                                # if(doRob && haveRob)

        varHat  <- rSqSum/(p*N)          # unbiased variance estimate - no Bessel correction
        corrFac <- 1/c4(p*N + 1)         # c4 correction with n = 2*N + 1
        chisqDF <- p*N                   # chi^2 degrees of freedom
    } else {                             # Singh C2 - true center is estimated
        rSqSum <- if(doRob && haveRob) { # sum squared radii centered data
            ## = N*trace of uncorrected covariance matrix -> N*trace((N-1)/N)*cov(xy)
            (N-1)*sum(diag(rob$cov))
        } else {
            xyCtr <- scale(xy, scale=FALSE, center=TRUE)  # centered data
            sum(xyCtr^2)                 # sum squared radii
        }                                # if(doRob && haveRob)

        ## unbiased variance estimate including Bessel correction
        varHat  <- rSqSum/(p*(N-1))
        corrFac <- 1/c4(p*N - (p-1))     # c4 correction with n = 2*N - 1
        chisqDF <- p*(N-1)               # chi^2 degrees of freedom
    }                                    # if(!missing(mu))

    sigHat  <- corrFac * sqrt(varHat)    # bias-corrected sigma estimate
    sigCIlo <- corrFac * sqrt(rSqSum / qchisq(1-(alpha/2), chisqDF))
    sigCIup <- corrFac * sqrt(rSqSum / qchisq(   alpha/2,  chisqDF))

    ## mean and sd for the chi distribution
    if(p == 1L) {                        # half normal with theta = sqrt(pi/2)
        ## http://mathworld.wolfram.com/Half-NormalDistribution.html
        MU <- sqrt(2/pi)
        SD <- sqrt((pi-2)/pi)
    } else if(p == 2L) {                 # Rayleigh distribution
        MU <- sqrt(pi/2)
        SD <- sqrt((4-pi)/2)
    } else if(p == 3L) {                 # Maxwell-Boltzmann distribution
        MU <- sqrt(8/pi)
        SD <- sqrt((3*pi-8)/pi)
    } else {
        MU <- NA_real_
        SD <- NA_real_
    }

    MR    <- sigHat * MU                 # radial error mean
    RSD   <- sigHat * SD                 # radial error standard deviation
    MRci  <- c(sigCIlo, sigCIup) * MU
    RSDci <- c(sigCIlo, sigCIup) * SD

    return(list(sigma=c(sigma=sigHat, sigCIlo=sigCIlo, sigCIup=sigCIup),
                  RSD=c(RSD=RSD, RSDciLo=RSDci[1], RSDciUp=RSDci[2]),
                   MR=c(MR=MR, MRciLo=MRci[1], MRciUp=MRci[2])))
}

#####---------------------------------------------------------------------------
## Rayleigh distribution
## pdf based on sigma
dRayleigh <-
function(x, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(x, scale)
    x     <- args[[1]]
    scale <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1L) { return(dens) }

    dens[keep] <- exp(log(x[keep]) - 0.5*(x[keep]/scale[keep])^2 - 2*log(scale[keep]))

    ## special case not caught so far
    dens[is.infinite(x)] <- 0
    return(dens)
}

## Rayleigh pdf based on sigma^2
dRayleighSq <-
function(x, scaleSq=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args    <- recycle(x, scaleSq)
    x       <- args[[1]]
    scaleSq <- args[[2]]

    dens <- numeric(length(x))
    keep <- which((x >= 0) | !is.finite(x))
    if(length(keep) < 1L) { return(dens) }

    dens[keep] <- exp(log(x[keep]) - 0.5*(x[keep]^2 / scaleSq[keep]) - log(scaleSq[keep]))

    ## special case not caught so far
    dens[is.infinite(x)] <- 0
    return(dens)
}

## Rayleigh cdf
pRayleigh <-
function(q, scale=1, lower.tail=TRUE) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(q, scale)
    q     <- args[[1]]
    scale <- args[[2]]

    pp   <- numeric(length(q))
    keep <- which((q >= 0) | !is.finite(q))

    if(lower.tail) {
        pp[keep] <- -expm1(-0.5 * (q[keep]/scale[keep])^2)
        ## some special values not caught before
        pp[which(q == -Inf)] <- 0
        pp[which(q ==  Inf)] <- 1
    } else {
        pp[keep] <- exp(-0.5 * (q[keep]/scale[keep])^2)
        ## some special values not caught before
        pp[which(q < 0)]    <- 1
        pp[which(q == Inf)] <- 0
    }

    return(pp)
}

## Rayleigh quantile function
qRayleigh <-
function(p, scale=1, lower.tail=TRUE) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)

    args  <- recycle(p, scale)
    p     <- args[[1]]
    scale <- args[[2]]

    keep <- which((p >= 0) & (p < 1))
    qq   <- rep(NA_real_, length(p))
    if(length(keep) < 1) { return(qq) }

    qq[keep] <- if(lower.tail) {
        scale[keep] * sqrt(-2 * log1p(-p[keep]))
    } else {
        scale[keep] * sqrt(-2 * log(p[keep]))
    }

    return(qq)
}

## random deviates
rRayleigh <-
function (n, scale=1) {
    is.na(scale) <- (scale <= 0) | !is.finite(scale)
    rn <- scale * sqrt(-2*log(runif(n)))
    return(rn)
}
