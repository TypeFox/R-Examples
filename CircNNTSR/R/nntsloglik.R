nntsloglik <-
function (data, cpars = 1/sqrt(2*pi), M = 0) 
{
    if (M == 0) 
        return(-length(data)*log(2 * pi))
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")
    y<-sum(log(nntsdensity(data, cpars,M)))
    return(y)
}

