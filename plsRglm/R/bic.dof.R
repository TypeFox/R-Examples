bic.dof <- function (RSS, n, DoF, sigmahat) 
{
    bic_temp <- RSS/n + log(n) * (DoF/n) * sigmahat^2
    return(bic_temp)
}
