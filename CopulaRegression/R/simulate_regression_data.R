
simulate_regression_data<-function(n,alpha,beta,R,S,delta,tau,family,zt=TRUE,exposure=rep(1,nrow(S))){
	mu<-as.vector(exp(R%*%alpha))
	lambda<-as.vector(exp(S%*%beta))*exposure
	#lambda<-E2lambda(m)
	mydata<-matrix(,n,2)
	for (i in 1:n){
    		mydata[i,]<-simulate_joint(1,mu[i],delta,lambda[i],theta=BiCopTau2Par(tau=tau,family=family),family,zt=zt)
	}
	return(mydata)
}

