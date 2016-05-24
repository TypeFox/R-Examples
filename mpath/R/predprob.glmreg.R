### derived from predprob.glm in pscl
predprob.glmreg <- function(obj, which, newdata = NULL, at = NULL, ...){
    newx <- newdata
    isNegBin <- obj$family=="negbin"
    isPoisson <- obj$family=="poisson"
    isBinomial <- obj$family=="binomial"
    if(!isNegBin & !isPoisson & !isBinomial)
        stop(paste("your object of class",class(obj),"is unsupported by predprob.glmreg"))
    if(is.null(newdata))
        yhat <- predict(obj,
                        type="response", which=which)
    else
        yhat <- predict(obj,
                        newx=newx,
                        type="response", which=which)

    y <- obj$y
    yUnique <- if(is.null(at)) 0:max(y) else at
    nUnique <- length(yUnique)
    p <- matrix(NA,length(yhat),nUnique)
    dimnames(p) <- list(NULL,yUnique)

    if(isNegBin){
        for(i in 1:nUnique){
            p[,i] <- dnbinom(mu=yhat,
                             size=obj$theta[which],
                             x=yUnique[i])
        }
    }

    if(isPoisson){
        for(i in 1:nUnique){
            p[,i] <- dpois(lambda=yhat,
                           x=yUnique[i])
        }
    }

    if(isBinomial){
        if(is.null(newdata))
            p <- predict(obj,
                            type="response", which)
        else
            p <- predict(obj,
                            newx=newx,
                            type="response", which)
        p <- cbind(1-p,p)
        dimnames(p) <- list(NULL,c("0","1"))
    }

    p
    
}

