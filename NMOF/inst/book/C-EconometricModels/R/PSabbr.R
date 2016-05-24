# PSabbr.R -- version 2010-12-30 
function (OF, algo = list(), ...) {
    # add defaults to algo, do a few sanity checks
    # ...
    mRU <- function(m,n) array(runif(m * n), dim = c(m, n))
    mRN <- function(m,n) array(rnorm(m * n), dim = c(m, n))
    d <- length(algo$max); vF <- numeric(algo$nP); vF[] <- NA
    mP <- algo$min + diag(algo$max-algo$min) %*% mRU(d,algo$nP)
    mV <- algo$initV * mRN(d,algo$nP)
    for (s in 1:algo$nP) vF[s] <- OF(mP[, s], ...)
    mPbest <- mP; vFbest <- vF
    sGbest <- min(vFbest); sgbest <- which.min(vFbest)[1]
    for (g in 1:algo$nG) {
        mDV <- algo$c1*mRU(d,algo$nP) * (mPbest - mP) + 
               algo$c2*mRU(d,algo$nP) * (mPbest[,sgbest] - mP)
        mV <- algo$iner * mV + mDV
        logik <- mV > 0
        mV[logik] <- pmin(mV, algo$maxV)[logik]
        logik <- mV < 0
        mV[logik] <- pmax(mV, -algo$maxV)[logik]
        mP <- mP + mV
        for (s in 1:algo$nP) vF[s] <- OF(mP[, s], ...)
        logik <- vF < vFbest
        mPbest[, logik] <- mP[, logik]
        vFbest[logik] <- vF[logik]
        if (min(vF) < sGbest) {
            sGbest <- min(vF); sgbest <- which.min(vF)[1]
        }
    }
    return(list(vPar = mPbest[,sgbest], 
                    OFvalue = sGbest, popF = vFbest))
}