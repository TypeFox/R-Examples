ztp.glm <-
function(y,S,exposure=rep(1,length(y)),sd.error=FALSE){
	foo<-function(beta,y,S,exposure){
		lambda<-exposure*as.vector(exp(S%*%beta))
		#out<-0
		#for (i in 1:length(y)){
		#	out<-out+log(dztp(y[i],lambda[i]))
		#}
		out<-sum(log(dztp(y,lambda)))
		return(-out)
	}
	beta0<- coef(glm(y ~ S-1, offset=log(exposure),family=poisson))
	#beta0<-c(log(mean(y)),rep(0,ncol(S)-1))
	sd=NULL
	my.model<-optim(beta0,foo,y=y,S=S,exposure=exposure,method="BFGS",hessian=sd.error)
	if (sd.error==TRUE){
		my.cov<-ginv(my.model$hessian)
		sd<-sqrt(diag(my.cov))
	}
	coefficients=my.model$par
	return(list(sd=sd,coefficients=coefficients))

}
