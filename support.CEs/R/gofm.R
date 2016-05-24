gofm <- function(output) 
{
# Name: gofm
# Title: Calculating goodness-of-fit measures
# Arguments
#  output   An object containing the output from the function clogit() or glm().


# set variables

    if (any(class(output) == "clogit") == TRUE) {     # clogit()
        N   <- output$nevent    # sample size
        LL0 <- output$loglik[1] # log likelihood value at start (all coefficiets are restrected to 0)
        LLb <- output$loglik[2] # log likelihood value at convergence
    } else if (any(class(output) == "glm") == TRUE) { # glm()
        N   <- nrow(output$data)
        LL0 <- -N * log(2)
        LLb <- as.vector(logLik(output))
    }
    K     <- length(output$coefficients) # number of estimated coefficients

# calculate various measures

    rho2  <- 1 - (LLb / LL0)              # rho-squared
    rho2a <- 1 - ((LLb - K) / LL0)        # adjusted rho-squared
    aic   <- -2 * LLb + 2 * K             # Akaike information criterion (AIC)
    bic   <- -2 * LLb + K * log(N)        # Bayesian information criterion (BIC)

# format and return output

    rtn <- list(RHO2    = rho2,
                AdjRHO2 = rho2a,
                AIC     = aic,
                BIC     = bic,
                K       = K,
                N       = N,
                LL0     = LL0,
                LLb     = LLb)

    class(rtn) <- "gofm"

    return(rtn)
}

