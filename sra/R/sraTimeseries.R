sraTimeseries <-
function (beta, delta = rep(0, length(beta)), mu0 = 0, logvarA0 = 0, 
    logvarE0 = 0, logNe = log(100), logn = log(1e+10), logvarM = log(1e-20), 
    kc = 0, kg = 0, o = mu0, s = 0) 
{
    ans <- list()
    if (is.na(o)) {
        o <- mu0
    }
    Ne <- exp(logNe)
    n <- exp(logn)
    varM <- exp(logvarM)
    mu <- mu0
    d <- 0
    vara <- exp(logvarA0)
    varA <- vara + d
    varE <- exp(logvarE0)
    for (t in 1:(length(beta))) {
        mu <- c(mu, mu[t] + varA[t] * (beta[t] + s * (o - mu[t])))
        deltaMu.t <- abs(mu[t] - o)
        deltaMu.tp1 <- abs(mu[t + 1] - o)
        vara.tp1 <- vara[t] * exp(kg * (deltaMu.tp1 - deltaMu.t)) * 
            (1 - 1/(2 * Ne)) + (1 - 1/Ne) * delta[t] * (varA[t]^2)/(2 * 
            n * (varA[t] + varE[t])) + varM
        d.tp1 <- 0.5 * (1 - 1/Ne) * (d[t] + (1 - 1/n) * delta[t] * 
            (varA[t]^2)/(varA[t] + varE[t]))
        varA.tp1 <- vara.tp1 + d.tp1
        varE.tp1 <- varE[t] * exp(kc * (deltaMu.tp1 - deltaMu.t))
        if (is.nan(varA.tp1) || is.infinite(varA.tp1)) {
            varA.tp1 <- 0
        }
        if (varA.tp1 < 0) {
            varA.tp1 <- 0
        }
        if (is.nan(varE.tp1) || is.infinite(varE.tp1)) {
            varE.tp1 <- 0
        }
        if (varE.tp1 < 0) {
            varE.tp1 <- 0
        }
        vara <- c(vara, vara.tp1)
        d <- c(d, d.tp1)
        varA <- c(varA, varA.tp1)
        varE <- c(varE, varE.tp1)
    }
    ans$mean <- mu
    ans$varA <- varA
    ans$varE <- varE
    ans$varP <- ans$varA + ans$varE
    return(ans)
}
