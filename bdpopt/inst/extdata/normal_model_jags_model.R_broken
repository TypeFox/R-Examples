## JAGS model file for the normal model

## All phase II and phase III responses are independent and normally distributed.
## There is one Emax model for the population mean of the efficacy response,
## and one Emax model for the population mean of the safety response.

## The Emax model used is
## Emax <- function(d, theta) {
##    theta[1] + theta[2] * d^theta[4] / (theta[3]^theta[4] + d^theta[4])
##}
## The typical parameter names maps to theta according to
## theta[1] = E0
## theta[2] = Emax
## theta[3] = ED50
## theta[4] = h
## d is dose

## The parameters for the efficacy Emax model are theta[1:4].
## The parameters for the safety Emax model are eta[1:4].

## k.II is the number of different doses in phase II.
## n.II is a vector of sample sizes of length k.II corresponding to the different doses in phase II.
## YE.II is the efficacy response vector of length k.II in phase II.
## YS.II is the safety response vector of length k.II in phase II.

## In phase III, there may be several independent trials, with YE.III[i] and YS.III[i] corresponding to the
## efficacy and safety response for the i:th trial, respectively. The number of trials is given by k.III.
## However, all trials share the same dose d.III and sample size n.III.

model {
    ## Phase II prior.
    ## Use normal prior for the E0 and Emax parameters.
    ## Use a lognormal prior for the ED50 and h parameters.
    ## In each case, mu is the mean parameter and tau the precision parameter.

    theta[1] ~ dnorm(theta.mu[1], theta.tau[1])
    theta[2] ~ dnorm(theta.mu[2], theta.tau[2])
    
    theta[3] ~ dlnorm(theta.mu[3], theta.tau[3])
    theta[4] ~ dlnorm(theta.mu[4], theta.tau[4])

    eta[1] ~ dnorm(eta.mu[1], eta.tau[1])
    eta[2] ~ dnorm(eta.mu[2], eta.tau[2])
    
    eta[3] ~ dlnorm(eta.mu[3], eta.tau[3])
    eta[4] ~ dlnorm(eta.mu[4], eta.tau[4])  
    
    ## Phase II model
    for (i in 1:k.II) {
        ## Efficacy
        muE.II[i] <- theta[1] + theta[2] * d.II[i]^theta[4] / (theta[3]^theta[4] + d.II[i]^theta[4])
        tauE.II[i] <- n.II[i] / sigmaE^2
        YE.II[i] ~ dnorm(muE.II[i], tauE.II[i])

        ## Safety
        muS.II[i] <- eta[1] + eta[2] * d.II[i]^eta[4] / (eta[3]^eta[4] + d.II[i]^eta[4])
        tauS.II[i] <- n.II[i] / sigmaS^2
        YS.II[i] ~ dnorm(muS.II[i], tauS.II[i])        
    }       

    ## Phase III model
    ## Efficacy
    muE.III <- theta[1] + theta[2] * d.III^theta[4] / (theta[3]^theta[4] + d.III^theta[4])
    tauE.III <- n.III / sigmaE^2    
    for (i in 1:k.III) {
        YE.III[i] ~ dnorm(muE.III, tauE.III)
    }

    ## Safety
    muS.III <- eta[1] + eta[2] * d.III^eta[4] / (eta[3]^eta[4] + d.III^eta[4])
    tauS.III <- n.III / sigmaS^2
    for (i in 1:k.III) {
        YS.III[i] ~ dnorm(muS.III, tauS.III)        
    }   
}
