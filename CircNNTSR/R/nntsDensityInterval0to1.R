nntsDensityInterval0to1 <-
function (S, cpars = 1/sqrt(2*pi), M = 0) 
{
    n<-length(S)
    if (M == 0) 
	return(rep(1,n))
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")
#   y <- 0
#   for (k in 0:M){
#        y <- y + cpars[k+1]*exp(1i * k * (2*pi) * S)
#   }
#   y <- y * Conj(y)
#   y <- (2*pi)*y
    y<-2*pi*nntsdensity(2*pi*S,cpars,M)
    return(y)
}

