Ffuncshifth<-function(time,t,lambda,mu) {
	i <- inter(time,t)
	if (time==0){
		res<-0
	} else {
	tnew<-c(t[1:i],time)
	res<-lambda[i]/(lambda[i]-mu[i])
	if (i>1) {
	for (j in 1:(i-1)) {
		res<-res*exp((lambda[j]-mu[j])*(tnew[j+1]-tnew[j]))}}
	res<-res*(exp((lambda[i]-mu[i])*(tnew[i+1]-tnew[i]))-1)
	res<-res+Ffuncshifth(t[i],t,lambda,mu)
	}
	res
}