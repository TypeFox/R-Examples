snntssimulation <- function(nsim=1,cpars=(1/(2*pi))^2,M=c(0,0))
{
    R <- 2
    if (R != length(M)) 
        return("Error: Dimensions of M and vector of observations are not equal")
    size <- length(cpars)
    if (size != prod(M + 1)) 
        return("Error: Length of cpars must be equal to prod(M+1)")
    if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")


    res<-matrix(0,nrow=nsim,ncol=2)
    conteo<-1
    for (k in 1:nsim){
	U1<-runif(1,0,2*pi)
	U2<-runif(1,0,pi)
	U3<-runif(1,0,prod(M+1)/((2*pi)^2))
	while (U3 > mnntsdensity(t(c(U1,U2)),cpars,M,R=2)){
		U1<-runif(1,0,2*pi)
		U2<-runif(1,0,pi)
		U3<-runif(1,0,prod(M+1)/((2*pi)^2))
		conteo<-conteo+1
	}
	res[k,]<-c(U1,U2)
    }
    resf<-list(simulations=res,conteo=conteo)
    resf
}







