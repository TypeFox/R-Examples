sraAutoregEvolvTimeseries <-
function (beta, delta = rep(0, length(beta)), mu0 = 0, logIA0 = 0, 
    logIE0 = 0, relativekA0 = 0, kA1 = 1, kA2 = 0, kA3 = 0, relativekE0 = 0, 
    kE1 = 1, kE2 = 0, kE3 = 0, threshold = 1e-10, logrelativekA0 = NULL, 
    logrelativekE0 = NULL, logkA1 = NULL, logkE1 = NULL, logkA2 = NULL, 
    logkE2 = NULL, logkA3 = NULL, logkE3 = NULL) 
{
    ans <- list()
    if (mu0 < threshold) {
        mu0 <- exp(mu0)/(exp(threshold))
    }
    mu <- mu0
    IA <- exp(logIA0)
    IE <- exp(logIE0)
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
    kA0 <- relativekA0 * IA
    kE0 <- relativekE0 * IE
    varA <- IA * mu * mu
    varE <- IE * mu * mu
    for (t in 1:(length(beta))) {
        mu <- c(mu, mu[t] + beta[t] * IA[t] * mu[t]^2)
        IA.tp1 <- kA0 + kA1 * IA[t]
        IE.tp1 <- kE0 + kE1 * IE[t]
        if (is.nan(IA.tp1) || is.infinite(IA.tp1)) {
            IA.tp1 <- threshold
        }
        if (IA.tp1 < threshold) {
            IA.tp1 <- threshold
        }
        if (is.nan(IE.tp1) || is.infinite(IE.tp1)) {
            IE.tp1 <- threshold
        }
        if (IE.tp1 < threshold) {
            IE.tp1 <- threshold
        }
        if (is.nan(mu[t + 1])) {
            mu[t + 1] <- threshold
        }
        if (is.infinite(mu[t + 1])) {
            mu[t + 1] <- 1/threshold
        }
        if (mu[t + 1] < threshold) {
            mu[t + 1] <- threshold
        }
        if (mu[t + 1] > 1/threshold) {
            mu[t + 1] <- (1/threshold)
        }
        if (t > 1) {
            IA.tp1 <- IA.tp1 + kA2 * IA[t - 1]
            IE.tp1 <- IE.tp1 + kE2 * IE[t - 1]
        }
        else {
            IA.tp1 <- IA.tp1 + kA2 * IA[1]
            IE.tp1 <- IE.tp1 + kE2 * IE[1]
        }
        if (t > 2) {
            IA.tp1 <- IA.tp1 + kA3 * IA[t - 2]
            IE.tp1 <- IE.tp1 + kE3 * IE[t - 2]
        }
        else {
            IA.tp1 <- IA.tp1 + kA3 * IA[1]
            IE.tp1 <- IE.tp1 + kE3 * IE[1]
        }
        IA <- c(IA, IA.tp1)
        IE <- c(IE, IE.tp1)
        varA <- c(varA, IA[t + 1] * (mu[t + 1]^2))
        varE <- c(varE, IE[t + 1] * (mu[t + 1]^2))
    }
    ans$mean <- mu
    ans$varA <- varA
    ans$varE <- varE
    ans$varP <- ans$varA + ans$varE
    return(ans)
}
