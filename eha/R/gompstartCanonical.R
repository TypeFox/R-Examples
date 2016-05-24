gompstartCanonical <- function(enter, exit, event, score){
    ## Gives start values for Gompertz parameters
    ## 
    ## _Only_ used in phreg (because of presence of 'score' here)
    ## This is 'canonical'

    ##
    ## This is for ONE stratum only! So input is from only one!

    ## Profiling; scale is profiled out:

    ## enter <- Y[, 1]
    ## exit <- Y[, 2]
    ## event <- Y[, 3]
    
    D <- sum(event)
    logD <- log(D)
    DT <- sum(exit * event)
    Dscore <- sum(event * score)
    
    shape.scale <- function(scale){
        escale <- exp(scale)
        eshape <- D / sum(exp(score) * (Hgompertz(exit, scale = escale, shape = 1,
                                    param = "canonical") -
                          Hgompertz(enter, scale = escale, shape = 1,
                                    param = "canonical")))
        shape <- log(eshape)
        shape
    }
        
    l.shape <- function(scale){
        n <- length(scale) # Vectorizing
        escale <- exp(scale)
        loglik <- numeric(n)
        for (i in 1:n){
            shape <- shape.scale(scale[i])
            loglik[i] <- D * (shape - scale[i]) + DT * exp(-scale[i]) + Dscore - D
        }
        loglik
    }

    scale <- log(max(exit)) # start value
    fit <- optim(scale, l.shape, control = list(fnscale = -1), method = "BFGS")
    scale <- fit$par
    shape <- shape.scale(scale)
    ret <- c(scale, shape)
    names(ret) <- c("scale", "shape")
    ret
}
    
