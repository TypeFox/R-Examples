.nll <- function(logLambda, y, m = 0, Ak, censorLimits)
{
    # Start with zero total
    nllTotal <- 0

    # Loop over antibodies
    for (i in names(y))
        nllTotal <- nllTotal + .nllGamma(logLambda = logLambda,
                                         y = y[, i],
                                         m = m,
                                         A = Ak$A[, i],
                                         k = Ak$k[, i],
                                         censorLimit = censorLimits[i])

    # Return total
    return(nllTotal)
}
