deltaGen <- function(y, X, offset = NULL, phiInit, thetaInit, type,
                     alpha, beta, alphaInit) {
    if (missing(beta)) {
        GLM <- initial(y, X, offset = offset, type = type, alpha = alpha)
        if (GLM$type == "NegBin") {
            delta <- c(GLM$beta, phiInit, thetaInit, GLM$alpha)
            names(delta) <- c(colnames(X), names(phiInit), names(thetaInit),
                              "alpha")
        } else {
            delta <- c(GLM$beta, phiInit, thetaInit)
            names(delta) <- c(colnames(X), names(phiInit), names(thetaInit))
        }
    } else {
        if (type == "NegBin") {
            if (length(beta) != ncol(X))
                stop("incorrect number of initial betas")
            if (missing(alphaInit))
                stop("initial alpha is not entered")
            delta <- c(beta, phiInit, thetaInit, alphaInit)
            names(delta) <- c(colnames(X), names(phiInit), names(thetaInit),
                              "alpha")
        } else {
            if (length(beta) != ncol(X))
                stop("incorrect number of initial betas")
            delta <- c(beta, phiInit, thetaInit)
            names(delta) <- c(colnames(X), names(phiInit), names(thetaInit))
        }
    }
    delta
}

thetaGen <- function(thetaLags, thetaInit) {
    if (missing(thetaLags) & !missing(thetaInit)) {
        stop("The specific orders are not entered")
    }
    if (missing(thetaLags) & missing(thetaInit)) {
        thetaLags <- numeric(0)
        thetaInit <- numeric(0)
    }
    if (!missing(thetaLags) & missing(thetaInit)) {
        thetaLags <- thetaLags
        thetaInit <- rep(0, length(thetaLags))
    }
    if (!missing(thetaLags) & !missing(thetaInit)) {
        if (length(thetaLags) != length(thetaInit))
            stop("the length of theta lags and theta initials differ")
        thetaLags <- thetaLags
        thetaInit <- thetaInit
    }
    names(thetaInit) <- paste(rep("theta", length(thetaInit)),
                              as.character(thetaLags), sep = "_")
    list(thetaLags = thetaLags, thetaInit = thetaInit)
}

phiGen <- function(phiLags, phiInit) {
    if (missing(phiLags) & !missing(phiInit)) {
        stop("The specific orders are not entered")
    }
    if (missing(phiLags) & missing(phiInit)) {
        phiLags <- numeric(0)
        phiInit <- numeric(0)
    }
    if (!missing(phiLags) & missing(phiInit)) {
        phiLags <- phiLags
        phiInit <- rep(0, length(phiLags))
    }
    if (!missing(phiLags) & !missing(phiInit)) {
        if (length(phiLags) != length(phiInit))
            stop("The length of theta lags and theta initials differ")
        phiLags <- phiLags
        phiInit <- phiInit
    }
    names(phiInit) <- paste(rep("phi", length(phiInit)), as.character(phiLags),
                            sep = "_")
    list(phiLags = phiLags, phiInit = phiInit)
}
