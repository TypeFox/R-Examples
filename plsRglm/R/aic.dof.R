aic.dof <- function (RSS, n, DoF, sigmahat) 
{
    aic_temp <- RSS/n + 2 * (DoF/n) * sigmahat^2
    return(aic_temp)
}
