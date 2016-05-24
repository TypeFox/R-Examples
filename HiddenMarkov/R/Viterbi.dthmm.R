Viterbi.dthmm <- function (object, ...){
    x <- object$x
    dfunc <- makedensity(object$distn)
    n <- length(x)
    m <- nrow(object$Pi)
    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    nu[1, ] <- log(object$delta) + dfunc(x=x[1], object$pm,
                                    getj(object$pn, 1), log=TRUE)
    logPi <- log(object$Pi)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi, 2, max) +
                      dfunc(x=x[i], object$pm, getj(object$pn, i),
                                  log=TRUE)
    }
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 
        1]] + nu[i, ])
    return(y)
}
