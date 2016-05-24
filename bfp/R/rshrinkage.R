#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[rshrinkage.R] by DSB Mon 05/10/2009 10:35 (CEST)>
##
## Description:
## Sample from the model-specific posterior of the shrinkage factor t = g / (1 + g).
##
## History:
## 04/09/2008   file creation
## 05/09/2008   correct innerFraction computation
## 05/10/2009   some beautifications
#####################################################################################


rshrinkage <- function(n=1,             # number of samples
                       R2,              # coefficient of determination in the model
                       nObs,            # number of observations
                       p,               # number of effects (without intercept)
                       alpha            # hyperparameter for hyper-g prior
                       )
{
    ## parameters for the beta functions
    shape1 <- (nObs - p - alpha + 1) / 2
    shape2 <- (p + alpha - 2) / 2

    ## uniforms needed for inverse sampling
    uniforms <- runif(n=n)
  
    ## evaluate the quantile function at these values to obtain the samples
    oneMinusR2 <- 1 - R2

    innerFraction <-
        oneMinusR2 /
            qbeta(uniforms +
                  (1 - uniforms) * pbeta(oneMinusR2, shape1, shape2),
                  shape1, shape2)
              
    ret <- (1 - innerFraction) / R2

    ## return the shrinkage samples
    return(ret)    
}
