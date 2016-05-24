"residuals.mmglm1" <-
function (object, ...) 
{
    m <- nrow(object$Pi)
    n <- length(object$y)
    tmp <- matrix(as.double(0), nrow=n, ncol=m)
    if (object$glmfamily$family=="binomial") size <- object$size
    else size <- rep(NA, n)
    #------------------------------------------
    #   calc forward-backward eqns
    for (k in 1:m)
        tmp[,k] <- dmmglm(object$y, object$beta[,k], object$sigma[k],
                          object$glmfamily, object$Xdesign, log=FALSE,
                          size=size)
    fb <- forwardback.dthmm(object$Pi, object$delta, tmp)
    #------------------------------------------
    #   calc conditional probability distribution
    for (i in 1:n)
        tmp[i,] <- pmmglm(object$y[i], object$beta, object$sigma, object$glmfamily,
                          object$Xdesign[i,], size=size[i])
    condprob <- probhmm(fb$logalpha, fb$logbeta, object$Pi, object$delta, tmp)
    #------------------------------------------
    #   discrete adjustment
    if (object$discrete==TRUE){
        for (i in 1:n)
            tmp[i,] <- pmmglm(object$y[i]-1, object$beta, object$sigma,
                              object$glmfamily, object$Xdesign[i,], size=size[i])
        condprob1 <- probhmm(fb$logalpha, fb$logbeta, object$Pi, object$delta, tmp)
        condprob <- (condprob + condprob1)/2
    }
    #------------------------------------------
    return(qnorm(condprob))
}

