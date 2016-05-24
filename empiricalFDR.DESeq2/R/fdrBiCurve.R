fdrBiCurve <-
function(fdr.table,maxLogP=5,...) {
	fdrt=fdr.table
	plot(dr~logp,fdrt,xlim=c(0.1,maxLogP),type="l",xaxt="n",mgp=c(2.3,1,0),xlab="pvalue",ylab="N(DEGs)",...)
	axis(1,labels=10^(-c(1:maxLogP)),at=c(1:maxLogP)) 
	lines(fdr~logp,fdrt,xlim=c(1,maxLogP),type="l",col="red")
	abline(v=-log(empiricalFDR(fdrt),10),lty=3)
	legend("topright",c("real","simul"),lty=1,col=c("black","red"),bty="n")
}
