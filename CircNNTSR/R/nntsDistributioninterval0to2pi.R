nntsDistributioninterval0to2pi <-
function (theta, cpars = 1/sqrt(2*pi), M = 0) 
{
    theta<-theta%%(2*pi)
    if (M == 0) 
        return(theta/(2*pi))
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")
#    y <- theta/(2*pi)
#    for (k in 0:M){
#	for (m in 0:M){
#		if (k != m)
#		y <- y + cpars[k+1]*Conj(cpars[m+1])*(1i/(k-m))*(1 - exp(1i*(k-m)*theta))
#		}
#   	 }
#    y <- Re(y)
    y<-nntsDistribution(theta,cpars,M)
    return(y)
}

