gmdl.dof <- function (sigmahat, n, DoF, yhat) 
{
    SS <- sigmahat^2
    denominator <- DoF * SS
    FF <- (yhat)/(DoF * SS)
    FF[FF == 0] = Inf
    gmdl_temp <- (n/2) * log(SS) + (DoF/2) * log(FF) + (1/2) * log(n)
    return(gmdl_temp)
}
