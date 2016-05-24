"Estep" <-
function (x, Pi, delta, distn, pm, pn = NULL) 
{
    dfunc <- makedensity(distn)
    m <- nrow(Pi)
    n <- length(x)
    y <- forwardback(x, Pi, delta, distn, pm, pn)
    logbeta <- y$logbeta
    logalpha <- y$logalpha
    LL <- y$LL
    u <- exp(logalpha + logbeta - LL)
    v <- array(NA, dim = c(n - 1, m, m))
    for (k in 1:m) {
        logprob <- dfunc(x=x[-1], getj(pm, k),
                         getj(pn, -1), log=TRUE)
        logPi <- matrix(log(Pi[, k]), byrow = TRUE, nrow = n - 
            1, ncol = m)
        logPbeta <- matrix(logprob + logbeta[-1, k],
            byrow = FALSE, nrow = n - 1, ncol = m)
        v[, , k] <- logPi + logalpha[-n, ] + logPbeta - LL
    }
    v <- exp(v)
    return(list(u = u, v = v, LL = LL))
}

