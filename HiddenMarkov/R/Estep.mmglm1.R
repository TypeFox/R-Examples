Estep.mmglm1 <- function(object, fortran=TRUE){
    #    x is a mmglm1 object
    m <- nrow(object$Pi)
    n <- length(object$y)
    prob <- matrix(as.double(0), nrow=n, ncol=m)
    for (k in 1:m)
        prob[,k] <- dmmglm(object$y, object$beta[,k], object$sigma[k],
                           object$glmfamily, object$Xdesign, log=FALSE,
                           size=object$size)
    y <- forwardback.dthmm(object$Pi, object$delta, prob, fortran=fortran)
    logbeta <- y$logbeta
    logalpha <- y$logalpha
    LL <- y$LL
    u <- exp(logalpha + logbeta - LL)
    v <- array(NA, dim=c(n-1, m, m))
    for (k in 1:m) {
        logprob <- log(prob[-1, k])
        logPi <- matrix(log(object$Pi[, k]), byrow=TRUE,
                        nrow=n-1, ncol=m)
        logPbeta <- matrix(logprob + logbeta[-1, k],
            byrow=FALSE, nrow=n-1, ncol=m)
        v[, , k] <- logPi + logalpha[-n,] + logPbeta - LL
    }
    v <- exp(v)
    return(list(u=u, v=v, LL=LL))
}

