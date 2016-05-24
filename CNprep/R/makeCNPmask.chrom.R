makeCNPmask.chrom<-function(imat,startcol=1,endcol=2,nprof=1,uthresh,dthresh){
	astart<-imat[,startcol]
	aend<-imat[,endcol]
	z<-cbind(c(astart,aend,aend+1),
		c(rep(1,length(astart)),rep(0,length(aend)),rep(-1,length(aend))))
	z<-z[order(z[,1]),]
	z[,2]<-cumsum(z[,2])
	z<-z[nrow(z)-rev(match(rev(unique(z[,1])),rev(z[,1])))+1,]
	#z[,1] gives unique start and end positions; z[,2] gives event counts there
	z<-cbind(z,z[,2]>=(uthresh*nprof)) #mark positions w/counts above upper thresh
	zsteps<-z[,3]-c(0,z[-nrow(z),3])
	ustart<-z[zsteps==1,1]
	zsteps<-z[,3]-c(z[-1,3],0)
	uend<-z[zsteps==1,1] #starts and ends of intervals w/count above upper thresh
	z[,3]<-z[,2]>=(dthresh*nprof)
	zsteps<-z[,3]-c(0,z[-nrow(z),3])
	dstart<-z[zsteps==1,1]
	zsteps<-z[,3]-c(z[-1,3],0)
	dend<-z[zsteps==1,1] #likewise for the lower thresh
	if(length(ustart)>0){
		ci<-containment.indicator(ustart,uend,dstart,dend)
		return(matrix(ncol=2,data=c( dstart[ci[,2]>=ci[,1]],dend[ci[,2]>=ci[,1]]),
			dimnames=list(NULL,c("start","end"))))
	} #ie intervals above lower thresh with counts above upper thresh inside
	else{
		return(matrix(ncol=3,nrow=0,dimnames=list(NULL,c("chrom","start","end"))))
	}
	}
