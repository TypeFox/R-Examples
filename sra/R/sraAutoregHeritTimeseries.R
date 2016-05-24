sraAutoregHeritTimeseries <-
function (beta, delta = rep(0, length(beta)), mu0 = 0, logith20 = 0, 
    logvarP0 = 0, relativekA0 = 0, kA1 = 1, kA2 = 0, kA3 = 0, 
    relativekE0 = 0, kE1 = 1, kE2 = 0, kE3 = 0, threshold = 1e-10, 
    logrelativekA0 = NULL, logrelativekE0 = NULL, logkA1 = NULL, 
    logkE1 = NULL, logkA2 = NULL, logkE2 = NULL, logkA3 = NULL, 
    logkE3 = NULL) 
{
    ans <- list()
    mu <- mu0
    h2 <- exp(logith20)/(1 + exp(logith20))
    varP <- exp(logvarP0)
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
    varA <- h2 * varP
    varE <- varP - varA
    h2 <- varA/(varA + varE)
    kA0 <- relativekA0 * h2
    kE0 <- relativekE0 * varP
    for (t in 1:(length(beta))) {
        mu <- c(mu, mu[t] + beta[t] * varA[t])
        h2.tp1 <- kA0 + kA1 * h2[t]
        varP.tp1 <- kE0 + kE1 * varP[t]
        threshold <- 1e-10
        if (is.nan(h2.tp1) || is.infinite(h2.tp1)) {
            h2.tp1 <- threshold
        }
        if (h2.tp1 < threshold) {
            h2.tp1 <- threshold
        }
        if (h2.tp1 > 1 - threshold) {
            h2.tp1 <- 1 - threshold
        }
        if (is.nan(varP.tp1) || is.infinite(varP.tp1)) {
            varP.tp1 <- 0
        }
        if (varP.tp1 < threshold) {
            varP.tp1 <- 0
        }
        if (t > 1) {
            h2.tp1 <- h2.tp1 + kA2 * h2[t - 1]
            varP.tp1 <- varP.tp1 + kE2 * varP[t - 1]
        }
        else {
            h2.tp1 <- h2.tp1 + kA2 * h2[1]
            varP.tp1 <- varP.tp1 + kE2 * varP[1]
        }
        if (t > 2) {
            h2.tp1 <- h2.tp1 + kA3 * h2[t - 2]
            varP.tp1 <- varP.tp1 + kE3 * varP[t - 2]
        }
        else {
            h2.tp1 <- h2.tp1 + kA3 * h2[1]
            varP.tp1 <- varP.tp1 + kE3 * varP[1]
        }
        h2 <- c(h2, h2.tp1)
        varP <- c(varP, varP.tp1)
        varA <- c(varA, h2.tp1 * varP.tp1)
        varE <- c(varE, varP.tp1 - varA[t + 1])
    }
    ans$mean <- mu
    ans$varA <- varA
    ans$varE <- varE
    ans$varP <- ans$varA + ans$varE
    return(ans)
}
