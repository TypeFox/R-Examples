sraEpiTimeseries <-
function (beta, delta = rep(0, length(beta)), mu0 = 0, logvarA0 = 0, 
    logvarE0 = 0, logNe = log(1000), logvarM = log(1e-20), logepsilon = 0, 
    logminusepsilon = NA, logvarepsilon = 0) 
{
    ans <- list()
    if (is.na(logepsilon) && is.na(logminusepsilon)) {
        stop("Either positive or negative epsilon must be considered")
    }
    if (!is.na(logepsilon)) {
        logminusepsilon <- NA
    }
    mu <- mu0
    varA <- exp(logvarA0)
    varE <- exp(logvarE0)
    varepsilon <- exp(logvarepsilon)
    if (!is.na(logepsilon)) {
        epsilon <- exp(logepsilon)
    }
    else {
        epsilon <- -exp(logminusepsilon)
    }
    varM <- exp(logvarM)
    Ne <- exp(logNe)
    for (t in 1:(length(beta))) {
        mu <- c(mu, mu[t] + beta[t] * (varA[t]))
        varA.tp1 <- varA[t] * (1 - 1/(2 * Ne)) + 2 * beta[t] * 
            epsilon * (varA[t]^2) + varM
        varE.tp1 <- varE[t] + (epsilon^2 + varepsilon) * (varA[t]^2)
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
        varA <- c(varA, varA.tp1)
        varE <- c(varE, varE.tp1)
    }
    ans$mean <- mu
    ans$varA <- varA
    ans$varE <- varE
    ans$varP <- ans$varA + ans$varE
    return(ans)
}
