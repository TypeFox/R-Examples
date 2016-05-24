fastGBM <-
function(Spot=1, sigma=0.2, n=1000, m=365, r=0.06, dr=0.0, mT=1) {
	GBM<-matrix(NA, nrow=n, ncol=m)
	for(i in 1:n) {
		GBM[i,]<-Spot*exp(cumsum(((r-dr)*(mT/m)-0.5*sigma*sigma*(mT/m))+(sigma*(sqrt(mT/m))*rnorm(m, mean= 0, sd=1))))
		}
	return(GBM)
}
