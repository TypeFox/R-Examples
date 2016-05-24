surfaceAICPlot<-function(fwd=NULL, bwd=NULL, out=NULL, summ=NULL, traitplot="none", cols=NULL, daic=FALSE, ...){

if(is.null(summ)){
	if(!is.null(out)){
		fwd<-out[[1]]
		bwd<-out[[2]]
	}
	summ<-surfaceSummary(fwd, bwd)
}

aics<-summ$aics
nreg<-summ$n_regimes_seq[2,]

if(traitplot%in%c("aic","dev")){
	nt<-dim(summ$lnls)[1]
	if(is.null(cols))cols<-rainbow(nt)
	lnls<-summ$lnls
	devs<--2*lnls
	if(traitplot=="aic"){
		#divide the "penalty component of AIC (from np & n) among traits
		resaic<-matrix(aics-colSums(devs),nrow=nt,ncol=length(aics),byrow=TRUE)/nt
		vals<-devs+resaic
	}else{
		vals<-devs
	}
	vals<-vals-vals[,1]
	aics<-aics-aics[1]
}else{
	vals<-NULL
	if(daic)aics<-aics-aics[1]
	}

	plot(nreg,aics,type="n",ylim=range(c(aics,vals)),...)
	abline(h=aics[1],lty="longdash")
	points(nreg,aics,type="l", ...)
	points(nreg,aics, pch=21, bg="black", ...)
if(traitplot%in%c("aic","dev")){
	for(i in 1:dim(vals)[1]){
		points(nreg,vals[i,],col=cols[i],type="l", ...)
		points(nreg,vals[i,],bg=cols[i],pch=21, ...)
		}
	}	
}	
