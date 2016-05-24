plotHR <-
function(dat,tmin,tmax,N,h=NULL,alpha=1/5,bound=Inf){
	a<-(N*tmin):(N*tmax); a<-a/N; z<-c()
	for (k in a){z<-c(z , HR(dat,k,h,alpha,bound))}
	plot(a,z,type="l",main=paste("Estimator of the hazard rate function"),xlab="Time",ylab="Hazard rate")
}
