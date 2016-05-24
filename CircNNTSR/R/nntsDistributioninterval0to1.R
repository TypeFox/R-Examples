nntsDistributioninterval0to1 <-
function (theta, cpars = 1/sqrt(2*pi), M = 0) 
{
    theta<-theta%%1
    if (M == 0) 
        return(theta) else {
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("sum of the squared norms of componentes greater than condition")
#    y <- theta
#    for (k in 0:M){
#for (m in 0:M){
#if (k != m)
#y <- y + cpars[k+1]*Conj(cpars[m+1])*(1i/(k-m))*(1 - exp(1i*(k-m)*(2*pi)*theta))
#}
#    }
    y<-nntsDistribution(2*pi*theta,cpars,M)
#    y <- Re(y)
    return(y)
    }
}

