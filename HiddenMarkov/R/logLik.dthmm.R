"logLik.dthmm" <-
function(object, fortran=TRUE, ...){
    m <- nrow(object$Pi)
    n <- length(object$x)
    dfunc <- makedensity(object$distn)
    prob <- matrix(as.double(0), nrow=n, ncol=m)
    for (k in 1:m)
        prob[,k] <- dfunc(x=object$x, getj(object$pm, k), object$pn, log=FALSE)
    y <- forwardback.dthmm(object$Pi, object$delta, prob, fortran,
                           fwd.only=TRUE)$LL
    return(y)
}


