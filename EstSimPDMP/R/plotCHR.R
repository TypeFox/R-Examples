plotCHR <-
function(dat,tmin,tmax,N){
	a<-(N*tmin):(N*tmax); a<-a/N; z<-c()
	for (k in a){z<-c(z , CHR( dat , k ))}
	plot(a,z,type="l",main=paste("Nelson-Aalen estimator of the cumulative hazard rate function"),xlab="Time",ylab="Cumulative hazard rate")
}
