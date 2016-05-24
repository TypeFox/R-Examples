plotProfile <- function(profile,shadow=TRUE,xlim,ylim,col.sea,col.bottom,xlab,ylab,...){
	
	if (ncol(profile)!=4) stop("The profile object should have 4 columns (Longitude, Latitude, Kilomtric distance from start of path and Depth)")
	if (missing(xlim)) xlim <- range(profile[,3])
	if (missing(ylim)) ylim <- range(c(0,profile[,4]))
	if (missing(col.bottom)) col.bottom <- rgb(198,184,151,maxColorValue=255)
	if (missing(col.sea)) col.sea <- rgb(130,180,212,maxColorValue=255)
	if (missing(xlab)) xlab <- "Distance from start of transect (km)"
	if (missing(ylab)) ylab <- "Depth (m)"

	xlim[1] <- max(xlim[1],min(profile[,3]))
	xlim[2] <- min(xlim[2],max(profile[,3]))

	profile <- profile[,-(1:2)]
	end <- which.min(abs(profile[,1]-xlim[2]))
	start <- which.min(abs(profile[,1]-xlim[1]))
	profile[end,1] <- xlim[2]
	profile[start,1] <- xlim[1]
	profile <- profile[start:end,]

	profile.poly <- rbind(c(xlim[1],ylim[1]), profile, c(xlim[2],ylim[1]))
	
	plot(profile, type="n", pch=19, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	polygon(matrix(c(xlim[1],0,xlim[1],ylim[1],xlim[2],ylim[1],xlim[2],0),byrow=TRUE,ncol=2),col=col.sea)
	abline(h=0,lwd=0.5, lty=1)
	if (shadow) lines(profile, col=rgb(110,110,110,.8,maxColorValue=255),lwd=4)
	polygon(profile.poly, col=col.bottom)
	polygon(matrix(c(xlim[1],ylim[1],xlim[1],ylim[1]-1000,xlim[2]+100,ylim[1]-1000,xlim[2]+100,ylim[1]),byrow=TRUE,ncol=2),col="white",border=NA)

	box()

}