"residuals.dthmm" <-
function (object, ...) 
{
    m <- nrow(object$Pi)
    n <- length(object$x)
    tmp <- matrix(as.double(0), nrow=n, ncol=m)
    #------------------------------------------
    #   calc forward-backward eqns
    dfunc <- makedensity(object$distn)
    for (k in 1:m)
        tmp[,k] <- dfunc(x=object$x, getj(object$pm, k), object$pn)
    fb <- forwardback.dthmm(object$Pi, object$delta, tmp)
    #------------------------------------------
    #   calc conditional probability distribution
    pfunc <- makedistn(object$distn)
    for (i in 1:n)
        tmp[i,] <- pfunc(object$x[i], object$pm, getj(object$pn, i))
    condprob <- probhmm(fb$logalpha, fb$logbeta, object$Pi, object$delta, tmp)
    #------------------------------------------
    #   discrete adjustment
    if (object$discrete==TRUE){
        for (i in 1:n)
            tmp[i,] <- pfunc(object$x[i]-1, object$pm, getj(object$pn, i))
        condprob1 <- probhmm(fb$logalpha, fb$logbeta, object$Pi, object$delta, tmp)
        condprob <- (condprob + condprob1)/2
    }
    #------------------------------------------
    return(qnorm(condprob))
}

