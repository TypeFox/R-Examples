plot.frequency <-
function(x, ...) {
### plotting frequency from frequency object

	freq <- x[[1]]
	for(i in 2:length(x)) freq <- cbind(freq, x[[i]])
	freq <- as.data.frame(freq)
	row.names(freq) <- attr(x, "row.names")
	names(freq) <- names(x)
		
	dim.data <- dim(freq)
	num.sectors <- dim.data[1] - 1
	num.classes <- dim.data[2] - 2
	if(num.classes>1) {
		freq.cum <- freq[1:num.sectors,3:dim.data[2]]
		for(i in 2:num.classes) freq.cum[,i] <- freq.cum[,i] + freq.cum[,i-1]
	} else {
		freq.cum <- data.frame(freq[1:num.sectors,2])
		num.classes <- 1
	}
	
	sectors <- seq(0, 360-360/num.sectors, by=360/num.sectors)
	sectors <- sectors+90
	sector.edges <- sectors*pi/180
	sector.width <- sector.edges[2] - sector.edges[1]
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col.set <- plot.param$col
	else {
		if(num.classes==1) col.set <- c("#4575B4")
		else if(num.classes==2) col.set <- c("#4575B4", "#91BFDB")
		else if(num.classes>2 && num.classes<=11) {
			if(requireNamespace("RColorBrewer", quietly=TRUE)) col.set <- rev(RColorBrewer::brewer.pal(num.classes, "RdYlBu"))
			else col.set <- rev(rainbow(num.classes, start=0.0, end=0.7))
		} else col.set <- rev(rainbow(num.classes, start=0.0, end=0.7))
	}
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex-0.2
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
	else cex.leg <- cex-0.2
	if(any(names(plot.param)=="width.leg")) width.leg <- plot.param$width.leg
	else width.leg <- 0.15
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="circles")) circles <- seq(plot.param$circles[1], plot.param$circles[2], by=plot.param$circles[3])
	else circles <- seq(5, 5*(trunc(ceiling(max(freq.cum, na.rm=TRUE))/5)+1), by=5)
	if(any(names(plot.param)=="fg")) fg <- plot.param$fg
	else fg <- FALSE
	if(any(names(plot.param)=="pos.axis")) pos.axis <- plot.param$pos.axis
	else pos.axis <- 60
	if(any(names(plot.param)=="sec.space")) sec.space <- 1-plot.param$sec.space
	else sec.space <- 0.8
	if(any(names(plot.param)=="col.circle")) col.circle <- plot.param$col.circle
	else col.circle <- "gray45"
	if(any(names(plot.param)=="col.cross")) col.cross <- plot.param$col.cross
	else col.cross <- "gray45"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "gray45"
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="col.leg")) col.leg <- plot.param$col.leg
	else col.leg <- "black"
	if(any(names(plot.param)=="lwd.circle")) lwd.circle <- plot.param$lwd.circle
	else lwd.circle <- 0.7
	if(any(names(plot.param)=="lwd.cross")) lwd.cross <- plot.param$lwd.cross
	else lwd.cross <- 0.7
	if(any(names(plot.param)=="lty.circle")) lty.circle <- plot.param$lty.circle
	else lty.circle <- 2
	if(any(names(plot.param)=="lty.cross")) lty.cross <- plot.param$lty.cross
	else lty.cross <- 1
	if(any(names(plot.param)=="title.leg")) title.leg <- plot.param$title.leg
	else title.leg <- "Wind speed\n[m/s]"
	if(any(names(plot.param)=="border.leg")) border.leg <- plot.param$border.leg
	else border.leg <- col.set[1:num.classes]
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="col.border")) col.border <- plot.param$col.border
	else col.border <- NULL
	if(any(names(plot.param)=="lwd.border")) lwd.border <- plot.param$lwd.border
	else lwd.border <- 0.5
	
	if(is.numeric(width.leg)) if(width.leg>1) stop("'width.leg' must be a numeric value between 0 and 1")
	if(num.classes>1 && width.leg!=0) lo <- layout(matrix(1:2, 1, 2), widths=c(1, width.leg))
	
	# plot
	plot.new()
	par(mar=c(1,1,1,1), las=1)
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1] > pin[2]) xlim <- (pin[1]/pin[2]) * xlim
	else ylim <- (pin[2]/pin[1]) * ylim
	plot.window(xlim, ylim, "", asp=1)
	
	if(!fg) {
		for(c in num.classes:1) {
			plot.data <- c(tail(rev(as.vector(freq.cum[,c])), n=1), head(rev(as.vector(freq.cum[,c])), n=-1))
			
			for (i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.data[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				polygon(xlist, ylist, col=col.set[c], border=col.set[c], lwd=0.01)
			}
		}
		if(!is.null(col.border)) {
			freq.cum[is.na(freq.cum)] <- 0
			plot.max <- apply(freq.cum, 1, max)
			plot.max <- c(tail(rev(plot.max), n=1), head(rev(plot.max), n=-1))
			for(i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.max[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				lines(xlist, ylist, col=col.border, lwd=lwd.border)
			}
		}
	}
	
	circle.pts <- seq(0, 2*pi, length.out=360)
	pos.axis <- pi/2 - pos.axis*pi/180
	
	for(i in 1:length(circles)) {
		rad <- 0.9 * circles[i]/tail(circles, 1)
		circle.x <- cos(circle.pts)*rad
		circle.y <- sin(circle.pts)*rad
		lines(circle.x, circle.y, lty=lty.circle, lwd=lwd.circle, col=col.circle)
		text(cos(pos.axis)*rad, sin(pos.axis)*rad, circles[i], cex=cex.axis, col=col.axis)
	}
	     
	lines(c(-0.92, 0.92), c(0, 0), lty=lty.cross, lwd=lwd.cross, col=col.cross)
	lines(c(0, 0), c(0.92, -0.92), lty=lty.cross, lwd=lwd.cross, col=col.cross)
	text(0, -0.90, "S", pos=1, cex=cex.lab, col=col.lab)
	text(-0.90, 0, "W", pos=2, cex=cex.lab, col=col.lab)
	text(0, 0.90, "N", pos=3, cex=cex.lab, col=col.lab)
	text(0.90, 0, "E", pos=4, cex=cex.lab, col=col.lab)
	text(1, -1, "%", pos=2, cex=cex.axis, col=col.axis)
	
	if(fg) {
		for(c in num.classes:1) {
			plot.data <- c(tail(rev(as.vector(freq.cum[,c])), n=1), head(rev(as.vector(freq.cum[,c])), n=-1))
			
			for (i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.data[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				polygon(xlist, ylist, col=col.set[c], border=col.set[c], lwd=0.01)
			}
		}
		if(!is.null(col.border)) {
			freq.cum[is.na(freq.cum)] <- 0
			plot.max <- apply(freq.cum, 1, max)
			plot.max <- c(tail(rev(plot.max), n=1), head(rev(plot.max), n=-1))
			for(i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.max[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				lines(xlist, ylist, col=col.border, lwd=lwd.border)
			}
		}
	}
	
	if(num.classes>1 && width.leg!=0) {
		par(mar=c(0,0,0,0))
		plot(0, type="n", axes=FALSE, xlab="", ylab="")
		legend("left", legend=names(freq.cum), title=title.leg, fill=col.set[1:num.classes], xjust=0, bty=bty.leg, border=border.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
	}
}
