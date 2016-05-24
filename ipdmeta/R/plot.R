plot_ipdlme <- function(x,y,...){
	
	# SIDE-BY-SIDE FOREST PLOT OF RANDOM EFFECTS IN IPD LME MODEL
	alpha <- x@ranef
	int <- alpha[,1]
	trt <- alpha[,2]
	
	se.int <- sapply(x@vcov.ranef,function(x) sqrt(x[1,1]))
	se.trt <- sapply(x@vcov.ranef,function(x) sqrt(x[2,2]))
	
	lower.int <- int-1.96*se.int
	upper.int <- int+1.96*se.int
	
	xlim <- range(c(lower.int,upper.int))
	
	par(mfrow=c(1,2))
	
	plot(y=1:nrow(alpha),x=int,axes=FALSE,xlim=xlim,xlab="Intercept",ylab="",pch=19,...)
	abline(v=0)
	
	axis(1)
	
	if(missing(y)) y <- row.names(alpha)
	axis(2,at=1:nrow(alpha),labels=y,las=1,cex.axis=.8)
	
	segments(x0=lower.int,x1=upper.int,y0=1:nrow(alpha),y1=1:nrow(alpha))

	lower.trt <- trt-1.96*se.trt
	upper.trt <- trt+1.96*se.trt

	xlim <- range(c(lower.trt,upper.trt))
	
	plot(y=1:nrow(alpha),x=trt,axes=FALSE,xlim=xlim,xlab="Treatment",ylab="",pch=19,...)
	abline(v=0)
	
	segments(x0=lower.trt,x1=upper.trt,y0=1:nrow(alpha),y1=1:nrow(alpha))
	axis(1)

}

setMethod("plot","ipdlme",plot_ipdlme)