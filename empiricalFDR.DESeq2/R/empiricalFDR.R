empiricalFDR <-
function(fdr.table,FDR=0.1,maxLogP=5,plot=FALSE,span=0.1,...) {
	if (plot==TRUE) { 
		plot(FDR~logp,fdr.table,xaxt="n",xlim=c(0.1,maxLogP),mgp=c(2.3,1,0),xlab="pvalue",ylab="FDR",...)
		axis(1,labels=10^(-c(1:maxLogP)),at=c(1:maxLogP)) 
	}
	loe=loess(FDR~logp, fdr.table,span=span)
	xs=seq(min(fdr.table$logp),max(fdr.table$logp),l=5000)
	fdr.lim=xs[which(predict(loe,newdata=data.frame(logp=xs))>FDR)][length(which(predict(loe,newdata=data.frame(logp=xs))>FDR))]
	if (plot==TRUE) { 
		lines(xs,predict(loe,newdata=data.frame(logp=xs)),...)
		abline(h=0.1,lty=3)
		abline(v=fdr.lim,lty=3)
	}
	return(10^(-fdr.lim)) 
}
