makeCNPmask<-function(imat,chromcol=1,startcol=2,endcol=3,nprof=1,uthresh,dthresh){
	CNPmask<-by(imat,INDICES=as.factor(imat[,chromcol]),FUN=makeCNPmask.chrom,
		startcol=startcol,endcol=endcol,nprof=nprof,uthresh=uthresh,
		dthresh=dthresh,simplify=T)
	myCNPmask<-matrix(ncol=2,byrow=T,data=unlist(lapply(CNPmask,t)))
	myCNPmask<-cbind(unlist(lapply(1:length(unique(imat[,chromcol])),
		FUN=function(x) rep(as.numeric(names(CNPmask)[x]),nrow(CNPmask[[x]])))),myCNPmask)
	dimnames(myCNPmask)[[2]]<-c("chrom","start","end")
	return(myCNPmask)
	}
