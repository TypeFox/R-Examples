mnntssimulation <- function(nsim=1,cpars=1/(2*pi),M=c(0,0),R=2)
{
    if (R != length(M)) 
        return("Error: Length of M and number of dimensions are not equal")
    size <- length(cpars)
    if (size != prod(M + 1)) 
        return("Error: Length of cpars must be equal to prod(M+1)")
    if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")
    res<-matrix(0,nrow=nsim,ncol=R)
    conteo<-1
    for (k in 1:nsim){
	U<-runif(R,0,2*pi)
	U2<-runif(1,0,prod(M+1)/((2*pi)^R))
	while (U2 > mnntsdensity(t(U),cpars,M,R)){
		U<-runif(R,0,2*pi)
		U2<-runif(1,0,prod(M+1)/((2*pi)^R))
		conteo<-conteo+1
	}
	res[k,]<-U
    }
    resf<-list(simulations=res,conteo=conteo)
    resf
}









