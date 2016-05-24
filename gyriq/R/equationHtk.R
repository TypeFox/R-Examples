equationHtk <- function(Htk, tiInd, covTerm)
## Evaluates at point 'Htk' the equation used to obtain the value of the 
## step function of the transformation model with censored data at the jump
## corresponding to the first observed failure time.
## 
## 'tiInd' is an indicator of the individuals forming the risk set and 'covTerm'
## is the scalar product of covariates with regression parameter estimates.
{
    rskVec <- exp(covTerm + Htk)
    return(sum(rskVec[tiInd]) - 1)
}