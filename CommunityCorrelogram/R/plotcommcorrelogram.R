plot.commcorrelogram<-function(x,y,alpha=0.05,...){
	x<-na.omit(x@community.correlogram)
	par(mfcol=c(2,1),mar=c(4,4,2,2))
	xy<-x[,c(1,3)]
	sym<-as.numeric(!(x[,4]<alpha | x[,4]>(1-alpha)))*2+19
	plot(xy,pch=sym,xlab='lag distance',type='o',...)
	abline(h=0,lty=1,lwd=0.5)
	par(mar=c(5,4,1,2))
	xz<-x[,c(1,4)]
	plot(xz,xlab='lag distance',
			ylab='significance',pch=19,ylim=c(0,1))
	abline(h=alpha,lty=2,col='red')
	legend('topleft',legend=c('empirical','alpha'),
		   pch=c(19,-1),lty=c(-1,1),col=c('black','red'),bg='white')
	
}
