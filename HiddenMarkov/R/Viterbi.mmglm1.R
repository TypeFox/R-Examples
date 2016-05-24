Viterbi.mmglm1 <- function (object, ...){
    x <- object$y
    n <- length(x)
    m <- nrow(object$Pi)
    nu <- matrix(NA, nrow=n, ncol=m)
    y <- rep(NA, n)
    if (object$glmfamily$family=="binomial") size <- object$size
    else size <- rep(NA, length(x))
    nu[1,] <- log(object$delta) + dmmglm(x[1], object$beta, object$sigma,
               object$glmfamily, object$Xdesign[1,], size=size[1], log=TRUE)
    logPi <- log(object$Pi)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i-1,], nrow=m, ncol=m)
        nu[i,] <- apply(matrixnu + logPi, 2, max) +
                   dmmglm(x[i], object$beta, object$sigma, object$glmfamily,
                          object$Xdesign[i,], size=size[i], log=TRUE)
    }
    if (any(nu[n,] == -Inf)) 
        stop("Problems With Underflow")
    y[n] <- which.max(nu[n,])
    for (i in seq(n-1, 1, -1)) y[i] <- which.max(logPi[,y[i+1]] + nu[i,])
    return(y)
}


