#'Make a three panel plot of a 3D plant
#'
#'@description Three plots of a 3D plant: views from the east, south and from above. This is
#'a lame way to plot the plant, as the stems are always plotted on top (whether
#'or not they are visible). It is available for quick plotting, and for
#'\code{\link{ypreport}}, as it does not require the \code{rgl} package.
#'
#'See \code{\link{plot.plant3d}} for more advanced, high quality, plotting of
#'plants.
#'
#'This function plots the plant from above, east and west views. Stems are also
#'plotted, as opposed to the standard plot of a projected plant (see
#'\code{\link{projectplant}}).
#'
#'@param plant An object of class 'plant3d' (see \code{\link{constructplant}}).
#'@param side Which side to plot (can specify more than 1).
#'@param stems If TRUE, plots the stem sections (always on top, lame).
#'@param autopar If TRUE, tries to guess how to split up the plotting device.
#'@note This function is called by \code{\link{ypreport}}.
#'@author Remko Duursma
#'@seealso \code{\link{plot.plant3d}},\code{\link{projectplant}}
#'@keywords misc
#'@examples
#'
#'
#'# Toona australis from above
#'viewplot(toona, "above")
#'
#'
#'
#'@export
viewplot <- function(plant, side=c("east","south","above"), stems=TRUE, autopar=TRUE){

	side <- tolower(side)
	
	if(plant$inputformat == "Q")stems <- FALSE
	
	addbranches <- function(i1, i2){
		x <- plant
		for(i in 1:length(x$stems)){

			l <- x$stems[[i]]
			segments(x0=l$xyz$from[i1], x1=l$xyz$to[i1],
				y0=l$xyz$from[i2], l$xyz$to[i2], col="brown", lwd=1)
		}
		for(i in 1:length(x$branches)){
			l <- x$branches[[i]]
			segments(x0=l$xyz$from[i1], x1=l$xyz$to[i1],
				y0=l$xyz$from[i2], l$xyz$to[i2], col="brown", lwd=1)
		}
	}

	if(autopar){
		if(length(side) == 3 && all(par()$mfrow==c(1,1)))
			par(mfrow=c(2,2), pty='s')
		if(length(side) == 2 && all(par()$mfrow==c(1,1)))
			par(mfrow=c(1,2), pty='s')
	}
	
	if("south" %in% side){
		pp <- projectplant(plant, altitude=0, azimuth=180)
		plot(pp, leaffill=TRUE,
			xlab="X",ylab="Z",main="View from South",
			ylim=c(0,pp$viewbound$maxy))
		if(stems)addbranches(1,3)
	}
	
	if("east" %in% side){
		pp <- projectplant(plant, altitude=0, azimuth=90)
		plot(pp, leaffill=TRUE,
			xlab="Y",ylab="Z",main="View from East",
			ylim=c(0,pp$viewbound$maxy))
		if(stems)addbranches(2,3)
	}
	
	if("above" %in% side){
		plot(projectplant(plant, altitude=90, azimuth=180), leaffill=TRUE,
			xlab="X",ylab="Y",main="View from above")
		if(stems)addbranches(1,2)
	}

}




