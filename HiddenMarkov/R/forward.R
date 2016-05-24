"forward" <-
function(x, Pi, delta, distn, pm, pn = NULL){
    #   forward probabilities alpha_ij
    m <- nrow(Pi)
    n <- length(x)
    phi <- matrix(delta, nrow=1)
    logalpha <- matrix(rep(NA, m*n), nrow=n)
    lscale <- 0
    dfunc <- makedensity(distn)
    for (i in 1:n){
        if (i > 1) phi <- phi %*% Pi
        phi <- phi %*% diag(dfunc(x[i], pm, getj(pn, i)))
        sumphi <- sum(phi)
        phi <- phi/sumphi
        lscale <- lscale + log(sumphi)
        logalpha[i,] <- log(phi) + lscale
    }
    return(logalpha)
}

