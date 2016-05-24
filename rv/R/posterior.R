## ========================================================================
## function  -  short description
## ========================================================================
## Name:        
## Description: 
## Parameters:  
## Required:    none
## History:     2004-06-  : 
##

posterior <- function(obj, ...) {
  UseMethod("posterior")
} 

## Old code:
##    for (sim in 1:nsim){
##      sigma[sim] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
##      beta[sim,] <- mvrnorm (1, beta.hat, V.beta*sigma[sim]^2)
##    }
##

## ========================================================================
## posterior.lm
## ========================================================================

posterior.lm <- function (obj, ...) {
  ## Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  ##
  summ <- summary(obj)
  sigma.hat <- summ$sigma
  beta.hat <- summ$coef[,1]
  V.beta <- summ$cov.unscaled
  k <- summ$df[1]
  n <- k + summ$df[2]
  sigma <- sigma.hat*sqrt((n-k)/rvchisq(1,n-k))
  beta.0 <- as.vector(beta.hat) + sigma * rvnorm(mean=0, var=V.beta)
  names(beta.0) <- names(beta.hat)
  list(beta=beta.0, sigma=sigma)
}

# ========================================================================
# posterior.glm
# ========================================================================

posterior.glm <- function(obj, ...) {
  ## Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  ##
  summ <- summary (obj, correlation=TRUE)
  beta.hat <- summ$coef[,1]
  sd.beta <- summ$coef[,2]
  corr.beta <- summ$corr
  k <- summ$df[1]
  n <- k + summ$df[2]
  V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
  ## dimnames(beta) <- list (NULL, names(beta.hat))
  beta <- rvnorm(mean=beta.hat, var=V.beta)
  list(beta=beta, sigma=V.beta)
}

