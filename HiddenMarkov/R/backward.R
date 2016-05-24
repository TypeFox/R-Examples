"backward" <-
function(x, Pi, distn, pm, pn = NULL){
    #   backward probabilities beta_ij
    m <- nrow(Pi)
    n <- length(x)
    logbeta <- matrix(rep(NA, m*n), nrow=n)
    logbeta[n,] <- 0
    phi <- matrix(rep(1/m, m), ncol=1)
    lscale <- log(m)
    dfunc <- makedensity(distn)
    for (i in seq(n-1, 1, -1)){
        phi <- Pi %*% diag(dfunc(x=x[i+1], pm,
                           getj(pn, i+1))) %*% phi
        logbeta[i,] <- log(phi) + lscale
        sumphi <- sum(phi)
        phi <- phi/sumphi
        lscale <- lscale + log(sumphi)
    }
    return(logbeta)
}

