sraAutoregTimeseries <-
function (beta, delta = rep(0, length(beta)), mu0 = 0, logvarA0 = 0, 
    logvarE0 = 0, relativekA0 = 0, kA1 = 1, kA2 = 0, kA3 = 0, 
    relativekE0 = 0, kE1 = 1, kE2 = 0, kE3 = 0, threshold = 1e-10, 
    logrelativekA0 = NULL, logrelativekE0 = NULL, logkA1 = NULL, 
    logkE1 = NULL, logkA2 = NULL, logkE2 = NULL, logkA3 = NULL, 
    logkE3 = NULL) 
{
    ans <- list()
    mu <- mu0
    varA <- exp(logvarA0)
    varE <- exp(logvarE0)
    if (!is.null(logrelativekA0)) {
        relativekA0 <- exp(logrelativekA0)
    }
    if (!is.null(logrelativekE0)) {
        relativekE0 <- exp(logrelativekE0)
    }
    if (!is.null(logkA1)) {
        kA1 <- exp(logkA1)
    }
    if (!is.null(logkE1)) {
        kE1 <- exp(logkE1)
    }
    if (!is.null(logkA2)) {
        kA2 <- exp(logkA2)
    }
    if (!is.null(logkE2)) {
        kE2 <- exp(logkE2)
    }
    if (!is.null(logkA3)) {
        kA3 <- exp(logkA3)
    }
    if (!is.null(logkE3)) {
        kE3 <- exp(logkE3)
    }
    kA0 <- relativekA0 * varA
    kE0 <- relativekE0 * varE
    for (t in 1:(length(beta))) {
        mu <- c(mu, mu[t] + beta[t] * varA[t])
        varA.tp1 <- kA0 + kA1 * varA[t]
        varE.tp1 <- kE0 + kE1 * varE[t]
        if (t > 1) {
            varA.tp1 <- varA.tp1 + kA2 * varA[t - 1]
            varE.tp1 <- varE.tp1 + kE2 * varE[t - 1]
        }
        else {
            varA.tp1 <- varA.tp1 + kA2 * varA[1]
            varE.tp1 <- varE.tp1 + kE2 * varE[1]
        }
        if (t > 2) {
            varA.tp1 <- varA.tp1 + kA3 * varA[t - 2]
            varE.tp1 <- varE.tp1 + kE3 * varE[t - 2]
        }
        else {
            varA.tp1 <- varA.tp1 + kA3 * varA[1]
            varE.tp1 <- varE.tp1 + kE3 * varE[1]
        }
        if (is.nan(varA.tp1) || is.infinite(varA.tp1)) {
            varA.tp1 <- 0
        }
        if (varA.tp1 < threshold) {
            varA.tp1 <- 0
        }
        if (is.nan(varE.tp1) || is.infinite(varE.tp1)) {
            varE.tp1 <- 0
        }
        if (varE.tp1 < threshold) {
            varE.tp1 <- 0
        }
        varA <- c(varA, varA.tp1)
        varE <- c(varE, varE.tp1)
    }
    ans$mean <- mu
    ans$varA <- varA
    ans$varE <- varE
    ans$varP <- ans$varA + ans$varE
    return(ans)
}
