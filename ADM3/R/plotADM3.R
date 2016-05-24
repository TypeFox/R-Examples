`plotADM3` <- function(l) {
	data=l$data; r=l$report;
	for(x in 1:length(r[,1])) {
		d=data[data[,1]==r[x,1],]
		st=r[x,2]-20000; so=r[x,3]+20000;
		s=max(abs(max(d[d[,2]>=st & d[,3]<=so,4])), abs(min(d[d[,2]>=st & d[,3]<=so,4])))
		plot(d[d[,2]>=st & d[,3]<=so,2 ], d[d[,2]>=st & d[,3]<=so,4], ylim=c(-s,s), pch=20, xlab="Position", ylab="LogRatio", main=paste("Feature", x, sep=""))
		abline(h=0, col="blue", lwd=2); abline(v=r[x,2]); abline(v=r[x,3])	
		segments(r[x,2], r[x,4], r[x,3], r[x,4], col="red", lwd=3)
		scan("")
	}
}
