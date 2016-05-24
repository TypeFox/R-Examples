Ffuncshift<-function(time,t,lambda,mu,rho=1) {
	1-rho+rho*(Ffuncshifth(time,t,lambda,mu)+1)
}