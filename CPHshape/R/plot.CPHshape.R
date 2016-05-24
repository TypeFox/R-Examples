plot.CPHshape <- function(x, cex=0.5, xlim=range(x$h.range), ylim=range(x$h.val), xlab="", ylab="hazard function", ...){
	
	ranges 	<- 	x$h.range
	val		<-	x$h.val
	type	<-	x$type
	mode	<-	x$mode	
	
	if(type=="increasing"){
		yy		<-	val
		xx		<-	ranges
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		plot(xx1, yy, pch=19, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
		segments(xx1,yy,xx2,yy, ...)
		abline(v=max(ranges), col=grey(0.5), lty=2)
		} 
	if(type=="decreasing"){
		yy		<-	val
		xx		<-	ranges
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		plot(xx2,yy, pch=19, cex=cex, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
		segments(xx1,yy,xx2,yy, ...)
		} 
	if(type=="unimodal"){
		plot(0, 0, type = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
		m0		<-	mode
		k0		<-	which(ranges==m0)[1]
		if(k0!=1){
		yy		<-	val[1:(k0-1)]
		xx		<-	ranges[1:k0]
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		points(xx1, yy, pch=19, cex=cex, ...)
		segments(xx1,yy,xx2,yy, ...)}
		if(k0==1){segments(0,0,m0,0, ...)}
		if(k0!=length(ranges)){
		yy		<-	val[(k0+1):length(val)]
		xx		<-	ranges[(k0+1):length(ranges)]
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		points(xx2, yy, pch=19, cex=cex, ...)
		segments(xx1,yy,xx2,yy, ...)}
		abline(v=m0, col=grey(0.5), lty=2)
		} 
	if(type=="ushaped"){
		plot(0, 0, type = 'n', xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
		a0		<-	mode
		k0		<-	max(which(ranges<a0))
		if(k0>=2){
		yy		<-	val[1:(k0-1)]
		xx		<-	ranges[1:k0]
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		points(xx1, yy, pch=19, cex=cex, ...)
		segments(xx1,yy,xx2,yy, ...)}
		if(k0==1){segments(0,0,ranges[2],0, ...)}
		if(k0<length(val)){
		yy		<-	val[(k0+1):length(val)]
		xx		<-	ranges[(k0+1):length(ranges)]
		xx1		<-	xx[-length(xx)]
		xx2		<-	xx[-1]
		points(xx2, yy, pch=19, cex=cex, ...)
		segments(xx1,yy,xx2,yy, ...)}
		if(k0==length(val)){segments(ranges[length(ranges)-1],0,ranges[length(ranges)],0, ...)}
		} 
	
	}	# plot.CPHshape

