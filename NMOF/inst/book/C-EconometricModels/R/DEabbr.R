# DEabbr.R -- version 2010-12-30
function (OF, algo = list(), ...) {
    # add defaults to algo, do a few sanity checks
    # ...
    mRU <- function(m, n) array(runif(m * n), dim = c(m, n))
    mRN <- function(m, n) array(rnorm(m * n), dim = c(m, n))
    shift <- function(x) c(x[length(x)], x[1:(length(x) - 1)])
    d  <- length(algo$max)
    vF <- numeric(algo$nP)
    vF[] <- NA; vFv <- vF
    mP <- algo$min + 
          diag(algo$max - algo$min) %*% mRU(d, algo$nP)
    for (s in 1:algo$nP) vF[s] <- OF(mP[ ,s], data = data)
    for (g in 1:algo$nG) {
        vI <- sample(1:algo$nP, algo$nP)
        R1 <- shift(vI); R2 <- shift(R1); R3 <- shift(R2)
        mPv <- mP[ ,R1] + algo$F * (mP[ ,R2] - mP[ ,R3])
        mI <- mRU(d, algo$nP) > algo$CR
        mPv[mI] <- mP[mI]
        for (s in 1:algo$nP) vFv[s] <- OF(mPv[ ,s], data = data)
        logik <- vFv < vF
        mP[ ,logik] <- mPv[ ,logik]
        vF[logik] <- vFv[logik]
        Fmat[g, ] <- vF
    }
    sGbest <- min(vF); sgbest <- which.min(vF)[1]
    return(list(vPar = mP[ ,sgbest], OFvalue = sGbest))
}
