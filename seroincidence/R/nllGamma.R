.nllGamma <- function(logLambda, y, m, A, k, censorLimit)
{
    n <- length(A)
    lambda <- exp(logLambda)
    log.rho <- ifelse(y > censorLimit,
                    sapply(X = y, FUN = .logRhoPdf, lambda = lambda, m = m, A = A, k = k, n = n), # uncensored
                    sapply(X = y, FUN = .logRhoCdf, lambda = lambda, m = m, A = A, k = k, n = n)) # censored

    return(-sum(log.rho))
}
