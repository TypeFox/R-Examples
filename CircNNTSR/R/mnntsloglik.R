mnntsloglik<-function (data, cpars = 1/sqrt(2 * pi), M = 0, R=1) 
{
    if (R != length(M)) 
        return("Error: Length of M and number of dimensions are not equal")
    size <- length(cpars)
    if (size != prod(M + 1)) 
        return("Error: Length of cpars must be equal to prod(M+1)")
    if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")
    if (sum(M) == 0) 
        return(-nrow(data) * R * log(2 * pi))
    y <- sum(log(mnntsdensity(data, cpars, M, R)))
    return(y)
}