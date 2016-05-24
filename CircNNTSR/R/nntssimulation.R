nntssimulation <- function(nsim=1,cpars=1/(2*pi),M=0)
{
    size <- length(cpars)
    if (size != M+1) 
        return("Length of cpars must be equal to M+1")
    if (abs(sum(Mod(cpars)^2) - 1/(2 * pi)) > 0.0000000001) 
        return("Sum of the squared norms of components greater than condition")
    res<-rep(0,nsim)
    conteo<-1
    for (k in 1:nsim){
	U1<-runif(1,0,2*pi)
	U2<-runif(1,0,(M+1)/(2*pi))
	while (U2 > nntsdensity(U1,cpars,M)){
		U1<-runif(1,0,2*pi)
		U2<-runif(1,0,(M+1)/(2*pi))
		conteo<-conteo+1
	}
	res[k]<-U1
     }
     resf<-list(simulations=res,conteo=conteo) 
     resf
}





