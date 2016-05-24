gompstartRate <- function(enter, exit, event, score){
    ## Gives start values for Gompertz parameters
    ## To be used only in phreg with 'rate' parametrization:
    ## h(t) = exp(shape + rate * t).
    ## score = X %*% init
    ##
    ## This is for ONE stratum only! So input is from only one!

    ## Profiling; shape is profiled out:

    ## enter <- Y[, 1]
    ## exit <- Y[, 2]
    ## event <- Y[, 3]
    
    D <- sum(event)
    logD <- log(D)
    DT <- sum(exit * event)

    shape.rate <- function(rate){
        eshape <- D * rate / sum(exp(score) * (exp(rate * exit) - exp(rate * enter)))
        log(eshape)
    }
        
    l.shape <- function(rate){
        ##cat("rate = ", rate, "\n")
        n <- length(rate) # Vectorizing
        loglik <- numeric(n)
        for (i in 1:n){
            shape <- shape.rate(rate[i])
            loglik[i] <- D * shape + DT * rate[i] + sum(event * score) -  D
        }
        loglik
    }

    rate <-  1 / max(exit) # start value
    fit <- optim(rate, l.shape, control = list(fnscale = -1), method = "BFGS")
    rate <- fit$par
    shape <- shape.rate(rate)
    ret <- c(rate, shape)
    names(ret) <- c("rate", "shape")
    ret
}
    
