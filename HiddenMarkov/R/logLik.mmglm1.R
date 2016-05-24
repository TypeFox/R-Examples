logLik.mmglm1 <- function(object, fortran=TRUE, ...){
    #    x is a mmglm1 object
    m <- nrow(object$Pi)
    n <- length(object$y)
    prob <- matrix(as.double(0), nrow=n, ncol=m)
    for (k in 1:m)
        prob[,k] <- dmmglm(object$y, object$beta[,k], object$sigma[k],
                           object$glmfamily, object$Xdesign, log=FALSE,
                           size=object$size)
    y <- forwardback.dthmm(object$Pi, object$delta, prob, fortran,
                           fwd.only=TRUE)$LL
    return(y)
}

