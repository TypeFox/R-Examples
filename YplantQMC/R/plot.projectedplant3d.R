#'@method plot projectedplant3d
#'@S3method plot projectedplant3d
#'@rdname projectplant
plot.projectedplant3d <- function(x, silhouette=FALSE, 
	xlim=NULL, ylim=NULL, leaffill=TRUE, leafcol="forestgreen",
	zerocenter=FALSE,xlab="X",ylab="Y",
	...){

	leaves <- lapply(x$leaves, function(m)data.frame(X=m$XYZ[,1], Y=m$XYZ[,2], Z=m$XYZ[,3]))
	
	# sort leaves by z value (if leaffill = TRUE, make sure 'closest' leaf is on top).
	zmean <- sapply(leaves, function(x)mean(x[,3]))
	leaves <- leaves[order(zmean)]
	
	xwid <-  x$viewbound$maxx -  x$viewbound$minx
	if(!is.null(xlim))xwid <- max(xwid, xlim[2] - xlim[1])
	ywid <-  x$viewbound$maxy -  x$viewbound$miny
	if(!is.null(ylim))ywid <- max(ywid, ylim[2] - ylim[1])
	span <- max(xwid, ywid)
    miny <- x$viewbound$miny
	minx <- x$viewbound$minx
	maxy <- x$viewbound$maxy
	maxx <- x$viewbound$maxx
	
	if(zerocenter)r <- max(abs(c(miny,minx,maxy,maxx)))
	
    if(is.null(ylim)){
		m <- mean(c(miny,maxy))
		ylim <- c(m-span/2, m+span/2)
	}
	if(is.null(xlim)){
		m <- mean(c(minx,maxx))
		xlim <- c(m-span/2, m+span/2)
	}
	
	if(!zerocenter)plot(1, type='n', xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab, ...)
	if(zerocenter){
		plot(1, type='n', xlim=c(-r,r),ylim=c(-r,r),xlab=xlab,ylab=ylab,...)
		abline(h=0,col="grey")
		abline(v=0,col="grey")
	}
	
	for(i in 1:length(leaves)){
		points(leaves[[i]], type='l')
		if(leaffill)polygon(x=leaves[[i]][,1], y=leaves[[i]][,2],col=leafcol)
	}
	
	if(silhouette){
	
		sil <- Silhouette(x)
		
		polygon(sil$xyz, border="grey")
		legend("topleft", paste("Silhouette area =",round(10^-2*sil$H,2),"cm2"),pch=-1,lty=-1,bty='n')
	}
	
}
