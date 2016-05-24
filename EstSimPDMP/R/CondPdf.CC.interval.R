CondPdf.CC.interval <-
function(dat,x,epsilon,tmin,tmax,nbre,h=NULL,alpha=1/5,verbose=TRUE,bound=Inf){
	
	tableau<-matrix(0,nrow=dim(dat)[1],ncol=2)
	tableau[,2]<-dat[,-1]
	for (i in 1:dim(tableau)[1]){
		if ( (dat[i,1]>x-epsilon) & (dat[i,1]<x+epsilon)){
			tableau[i,1]<-1
		} else {tableau[i,1]<-0}
	}
	print(paste("Effective sample size =",sum(tableau[,1])))
	
	z<-CondPdf.DC.interval(tableau,1,tmin,tmax,nbre,h,alpha,verbose=FALSE,bound)
	grid<-z$pdf
	a<-z$times
	
	if (verbose){
		if (length(x)>1){	
			plot(a,grid,type="l",main=paste("Estimator of the conditional density"),xlab="Time",ylab=paste("Conditional density given state=(",paste(x,collapse=","),")",sep=""))
			} else {
				plot(a,grid,type="l",main=paste("Estimator of the conditional density"),xlab="Time",ylab=paste("Conditional density given state=",x,sep=""))
			}
		}
	list(times=a,pdf=grid)
}
