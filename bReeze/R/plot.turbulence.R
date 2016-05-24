plot.turbulence <-
function(x, ...) {
### plotting turbulence intensity from turbulence object
	
	turb <- x[[1]]
	for(i in 2:length(x)) turb <- cbind(turb, x[[i]])
	turb <- as.data.frame(turb)
	row.names(turb) <- attr(x, "row.names")
	names(turb) <- names(x)
	
	num.sectors <- dim(turb)[1] - 1
	sectors <- seq(0, 360-360/num.sectors, by=360/num.sectors)
	sectors <- sectors+90
	sector.edges <- sectors*pi/180
	sector.width <- sector.edges[2] - sector.edges[1]
	
	turb.max <- 100*max(turb$total[1:num.sectors], na.rm=TRUE)
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else col <- "#E41A1C"
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex-0.2
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="circles")) circles <- seq(plot.param$circles[1], plot.param$circles[2], by=plot.param$circles[3])
	else circles <- seq(5, 5*(trunc(ceiling(turb.max)/5)+1), by=5)/100
	if(any(names(plot.param)=="fg")) fg <- plot.param$fg
	else fg <- FALSE
	if(any(names(plot.param)=="col.circle")) col.circle <- plot.param$col.circle
	else col.circle <- "gray45"
	if(any(names(plot.param)=="col.cross")) col.cross <- plot.param$col.cross
	else col.cross <- "gray45"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "gray45"
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="lwd.circle")) lwd.circle <- plot.param$lwd.circle
	else lwd.circle <- 0.7
	if(any(names(plot.param)=="lwd.cross")) lwd.cross <- plot.param$lwd.cross
	else lwd.cross <- 0.7
	if(any(names(plot.param)=="lty.circle")) lty.circle <- plot.param$lty.circle
	else lty.circle <- 2
	if(any(names(plot.param)=="lty.cross")) lty.cross <- plot.param$lty.cross
	else lty.cross <- 1
	if(any(names(plot.param)=="pos.axis")) pos.axis <- plot.param$pos.axis
	else pos.axis <- 60
	if(any(names(plot.param)=="sec.space")) sec.space <- 1-plot.param$sec.space
	else sec.space <- 0.8
	if(any(names(plot.param)=="col.border")) col.border <- plot.param$col.border
	else col.border <- col
	if(any(names(plot.param)=="lwd.border")) lwd.border <- plot.param$lwd.border
	else lwd.border <- 0.5
	
	# plot
	par(mar=c(1,1,1,1), las=1)
	plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1] > pin[2]) xlim <- (pin[1]/pin[2]) * xlim
	else ylim <- (pin[2]/pin[1]) * ylim
	plot.window(xlim, ylim, "", asp=1)
	
	if(!fg) {
		plot.data <- c(tail(rev(as.vector(turb$total[1:num.sectors])), n=1), head(rev(as.vector(turb$total[1:num.sectors])), n=-1))    	
		for (i in 1:num.sectors) {
			arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
			rad <- 0.95 * plot.data[i]/tail(circles, 1)
			xlist <- c(0, rad * cos(arc.pts), 0)
			ylist <- c(0, rad * sin(arc.pts), 0)
		   	polygon(xlist, ylist, col=col, border=col.border, lwd=lwd.border)
		}
	}
	
   	circle.pts <- seq(0, 2*pi, length.out=360)
   	pos.axis <- pi/2 - pos.axis*pi/180
	
	for(i in 1:length(circles)) {
		rad <- 0.95 * circles[i]/tail(circles, 1)
		circle.x <- cos(circle.pts)*rad
		circle.y <- sin(circle.pts)*rad
		lines(circle.x, circle.y, lwd=lwd.circle, lty=lty.circle, col=col.circle)
		text(cos(pos.axis)*rad, sin(pos.axis)*rad, circles[i], cex=cex.axis, col=col.axis)
	}
	    
	lines(c(-0.97, 0.97), c(0, 0), lwd=lwd.cross, lty=lty.cross, col=col.cross)
	lines(c(0, 0), c(0.97, -0.97), lwd=lwd.cross, lty=lty.cross, col=col.cross)
	text(0, -0.95, "S", pos=1, cex=cex.lab, col=col.lab)
	text(-0.95, 0, "W", pos=2, cex=cex.lab, col=col.lab)
	text(0, 0.95, "N", pos=3, cex=cex.lab, col=col.lab)
	text(0.95, 0, "E", pos=4, cex=cex.lab, col=col.lab)
	
	if(fg) {
		plot.data <- c(tail(rev(as.vector(turb$total[1:num.sectors])), n=1), head(rev(as.vector(turb$total[1:num.sectors])), n=-1))    	
		for (i in 1:num.sectors) {
			arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
			rad <- 0.95 * plot.data[i]/tail(circles, 1)
			xlist <- c(0, rad * cos(arc.pts), 0)
			ylist <- c(0, rad * sin(arc.pts), 0)
		   	polygon(xlist, ylist, col=col, border=col.border, lwd=lwd.border)
		}
	}
}
