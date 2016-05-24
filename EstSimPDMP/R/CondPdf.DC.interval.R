CondPdf.DC.interval <-
function( dat , x , tmin , tmax , nbre , h=NULL , alpha=1/5 , verbose=TRUE , bound=Inf){
	
	N<-length(dat[1,])-1; m<-length(dat[,1]);
	
	A<-dat[,1:N]; B<-as.matrix(dat[ !duplicated(A) ,1:N]);
	a<-(nbre*tmin):(nbre*tmax); a<-a/nbre; z<-c()
	
	for (t in a){
		s<-0
		for (k in 1:length(B[,1])){
			y<-c();
			for (l in 1:N){
				y<-c(y , B[k,l])
			}
			D<-.Tri2( dat , x , y)
			if (length(D)>0){
				s<-s+HR(D , t , h , alpha , bound)*.CondSurv(dat , x , y , t)
				}
			}
		z=c(z,s)
		}
	if (verbose){
		if (length(x)>1){	
			plot(a,z,type="l",main=paste("Estimator of the conditional density"),xlab="Time",ylab=paste("Conditional density given state=(",paste(x,collapse=","),")",sep=""))
			} else {
				plot(a,z,type="l",main=paste("Estimator of the conditional density"),xlab="Time",ylab=paste("Conditional density given state=",x,sep=""))
			}
		}	
	list(times=a,pdf=z)	
}
