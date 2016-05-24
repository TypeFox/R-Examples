RUVrinv <-
function (Y, X, ctl, Z = 1, eta = NULL, fullW0 = NULL, invsvd = NULL, 
    lambda = NULL, k = NULL, l = NULL, randomization = FALSE, 
    iterN = 1e+05, inputcheck = TRUE) 
{
    if (inputcheck) 
        inputcheck1(Y, X, Z, ctl)
    Y = RUV1(Y, eta, ctl)
    if (!is.null(lambda)) 
        return(RUVinv(Y, X, ctl, Z = Z, fullW0 = fullW0, invsvd = invsvd, 
            lambda = lambda, randomization = randomization, iterN = iterN))
    if (is.null(k)) {
        if (is.null(l) & (ncol(X) > 1)) 
            warning("Neither lambda nor k are specified, so a call to getK will be made.  But p > 1 and l is not specified.  Arbitrarily setting l = 1.")
        temp = getK(Y, X, ctl, Z = Z, fullW0 = fullW0, inputcheck = FALSE)
        k = temp$k
        fullW0 = temp$fullW0
    }
    ruv4fit = RUV4(Y, X, ctl, k, Z = Z, fullW0 = fullW0)
    lambda = sum(ruv4fit$sigma2[ctl])
    return(RUVinv(Y, X, ctl, Z = Z, fullW0 = fullW0, invsvd = invsvd, 
        lambda = lambda, randomization = randomization, iterN = iterN))
}
