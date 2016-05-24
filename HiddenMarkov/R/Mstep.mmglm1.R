Mstep.mmglm1 <- function(object, u){
    m <- ncol(u)
    beta <- NULL
    sigma <- NULL
    if (object$glmfamily$family=="binomial")
        y <- cbind(object$y, object$size-object$y)
    else y <- object$y
    for (k in 1:m){
        w <- u[,k]
        z <- glm.fit(object$Xdesign, y, weights=w,
                     start=object$beta[,k], family=object$glmfamily,
                     control=glm.control(maxit=500))
        fitval <- fitted(z)
        var <- object$glmfamily$variance(fitval)
        if (object$glmfamily$family=="binomial"){
            fitval <- object$size*fitval
            var <- object$size*var
        }

#       if (distn!="binomial") sigma[k] <- sqrt(sum(w*(x-fitval)^2/var)/sum(w))
#       else sigma[k] <- sqrt(sum(w*(x[,1]/pn$size-fitval)^2/var)/sum(w))

        sigma <- c(sigma, sqrt(sum(w*(object$y-fitval)^2/var)/sum(w)))

#         if (object$glmfamily$family!="binomial")
#             sigma <- c(sigma, sqrt(sum(w*(object$y-fitval)^2/var)/sum(w)))
#         else sigma <- c(sigma, sqrt(sum(w*object$size*(object$y/object$size-fitval)^2/var)/sum(w)))

        beta <- cbind(beta, as.vector(coefficients(z)))
    }
    return(list(beta=beta, sigma=sigma))
}


