gs5 <-
function (sigma, x, y, distr, asigma, bsigma, delta) 
{
    sigmaStar <- rgamma(1, shape = delta, rate = delta/sigma)
    sigmaT <- sigma
    qgammas <- sigmaT/sigmaStar
    Qgammas <- sigmaStar/sigmaT
    Term2 <- qgammas^(2 * delta - 1) * exp(-delta * (qgammas - 
        Qgammas))
    Kgamma <- Qgammas^(asigma - 1) * exp(-bsigma * (sigmaStar - 
        sigmaT))
    Prod <- 1
    for (i in seq(length(x))) {
        Prod <- Prod * (dk(x[i], distr = distr, mu = y[i], sigma = sigmaStar)/dk(x[i], 
            distr = distr, mu = y[i], sigma = sigmaT))
    }
    q3 <- min(1, Kgamma * Prod * Term2)
    sigma <- ifelse(runif(1) <= q3, sigmaStar, sigmaT)
    return(sigma)
}
