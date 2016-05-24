IyenGreenMLE <- function(t, q, N, type = 1, alpha = 0.05){

    ## Compute MLE for weight function w_2 in Iyengar & Greenhouse (1988, p. 112)
    ## note that the parametrization here is that x lives on the 
    ## t-statistic scale (not on p-value scale)

    ## compute maximum likelihood estimates
    para0 <- c(1, 1)
    res1 <- optim(par = para0, fn = IyenGreenLoglikT, method = "BFGS", control = list(trace = TRUE, REPORT = 1, fnscale = -1), t = t, q = q, N = N, type = type)
    theta <- res1$par[1]
    beta <- res1$par[2]
        
    ## p-values
    p0 <- 2 * pnorm(-abs(t))    
        
    ## generate output object
    res <- list("theta" = theta, "beta" = beta, "p" = p0)    
    return(res)
    }

