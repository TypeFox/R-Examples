gompstartAft <- function(enter, exit, event){
    ## Gives start values for Gompertz parameters
    ## To be used in aftreg.
    ##
    ## This is for ONE stratum only! So input is from only one!

    ## NOTE: Revised in 2.3-2: Using the 'canonical' representation:
    ## h(t, (gamma, alpha)) = exp(gamma - alpha) * exp(t * exp(-alpha))
    
    ## Profiling; gamma is profiled out:
    D <- sum(event)
    logD <- log(D)
    DT <- sum(exit * event)
    funk <- function(alpha){
        n <- length(alpha) # Vectorizing
        loglik <- numeric(n)
        for (i in 1:n){
            ealpha <- exp(-alpha[i])
            S <- sum(exp(exit * ealpha) - exp(enter * ealpha))
        ##    gamma <- logD - alpha[i] - log(S)
            gamma <- logD - log(S) # 2.3-2
        ##    loglik[i] <- D * gamma + DT * ealpha - D
            loglik[i] <- D * (gamma - alpha[i]) + DT * ealpha - D # 2.3-2
        }
        loglik
    }

    alpha <- log(max(exit)) # start value
    ##from <- alpha - width / 2
    ##to <- alpha + width / 2
    fit <- optim(alpha, funk, control = list(fnscale = -1), method = "BFGS")
    alpha <- fit$par
    S <- sum(exp(exit * exp(-alpha)) - exp(enter * exp(-alpha)))
    ## gamma <- log(D) - alpha - log(S)
    gamma <- log(D) - log(S) # 2.3-2    
    ret <- c(alpha, gamma)
    ##names(ret) <- c("alpha", "gamma")
    ret
}
    
