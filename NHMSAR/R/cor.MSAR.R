cor.MSAR <-
function(data,data.sim,lag=NULL,nc=1,alpha=.05,plot=FALSE,xlab="Time (days)") {
	N.s = dim(data)[2]
	N.sim = dim(data.sim)[2]
	N = floor(N.sim/N.s)
	if (is.null(lag)) {lag = floor(dim(data)[1]/2)}
	C.data=matrix(0,lag,1)
	for (ex in 1:dim(data)[2]) {
		C = acf(data[,ex,nc], lag.max=lag-1,type = "correlation",plot=FALSE)
		C.data=C.data+c(C$acf)
	}
	C.data = C.data/N.s
	C.sim=matrix(0,lag,N)
	for (k in 1:N) {
		for (ex in ((k-1)*N.s+1):(k*N.s)) {
			C = acf(data.sim[,ex,nc], lag.max=lag-1, type = "correlation",plot=FALSE)
			C.sim[,k]=C.sim[,k]+c(C$acf)
		}
		C.sim[,k] = C.sim[,k]/N.s
	}
	IC = matrix(NA,2,lag)
	for (l in 1:lag) {
		IC[,l] = quantile(C.sim[l,],probs=c(alpha/2,1-alpha/2))
	}
	C.sim = apply(C.sim,1,mean)
	if (plot) {
		plot(0:(lag-1),C.data,typ="l",ylab="Correlation",xlab=xlab,lwd=2)
		lines(0:(lag-1),C.sim,col="red")
		lines(0:(lag-1),IC[1,],col="red",lty=3)
		lines(0:(lag-1),IC[2,],col="red",lty=3)
	}

	return(list(C.data=C.data,C.sim=C.sim,CI.sim=IC,lags=0:(lag-1)))
}
