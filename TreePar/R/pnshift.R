pnshift<-function(n,time,t,lambda,mu,rho=1){
	i <- inter(time,t)
	rho1<-rho
	rho<-lambda*0+1
	rho[1]<-rho1
	probext<-q2(i,time,t,lambda,mu,rho)
	if (n==0){res<-probext} else {
	res<-1-probext
	Finv<- 1/Ffuncshift(time,t,lambda,mu,rho1)
	res<-res*Finv*(1-Finv)^(n-1)}
	res
}