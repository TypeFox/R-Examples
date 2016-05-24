plot.energy <-
function(x, show.total=TRUE, ...) {
###	plotting wind energy rose
	
	en <- x[[1]]
	for(i in 2:length(x)) en <- cbind(en, x[[i]])
	en <- as.data.frame(en)
	row.names(en) <- attr(x, "row.names")
	names(en) <- names(x)
	
	dim.data <- dim(en)
	num.sectors <- dim.data[1] - 1
	num.classes <- dim.data[2] - 1
	total <- en$total[num.sectors+1]
	unit <- attr(x, "unit")
	
	if(num.classes>1) {
		e.cum <- en[1:num.sectors,2:dim.data[2]]
		for(i in 2:num.classes) e.cum[,i] <- e.cum[,i] + e.cum[,i-1]
	} else {
		e.cum <- data.frame(en[1:num.sectors,1])
		num.classes <- 1
	}
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	sectors <- seq(0, 360-360/num.sectors, by=360/num.sectors)
	sectors <- sectors+90
	sector.edges <- sectors*pi/180
	sector.width <- sector.edges[2] - sector.edges[1]
		
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else {
		if(num.classes==1) col <- c("#31A354")
		else if(num.classes==2) col <- c("#31A354", "#A1D99B")
		else if(num.classes>2 && num.classes<=11) {
			if(requireNamespace("RColorBrewer", quietly=TRUE)) col <- rev(RColorBrewer::brewer.pal(num.classes, "Greens"))
			else rev(rainbow(num.classes, start=0.0, end=0.7))
		} else col <- rev(rainbow(num.classes, start=0.0, end=0.7))
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
	else width.leg <- 0.2
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="circles")) circles <- seq(plot.param$circles[1], plot.param$circles[2], by=plot.param$circles[3])
	else circles <- NULL
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
	else border.leg <- col[1:num.classes]
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="col.border")) col.border <- plot.param$col.border
	else col.border <- NULL
	if(any(names(plot.param)=="lwd.border")) lwd.border <- plot.param$lwd.border
	else lwd.border <- 0.5
	
	if(is.null(circles)) {
		energy.max <- max(e.cum, na.rm=TRUE)
		mag <-length(strsplit(as.character(energy.max),"")[[1]])-1
		circ.max <- ceiling(energy.max/10^mag)*10^mag
		if(circ.max<=2*10^mag) circles <- tail(seq(0, 2*10^mag, length.out=5), -1)
		if(circ.max>2*10^mag && circ.max<=3*10^mag) circles <- tail(seq(0, 3*10^mag, length.out=4), -1)
		if(circ.max>3*10^mag && circ.max<=4*10^mag) circles <- tail(seq(0, 4*10^mag, length.out=5), -1)
		if(circ.max>4*10^mag && circ.max<=5*10^mag) circles <- tail(seq(0, 5*10^mag, length.out=5), -1)
		if(circ.max>5*10^mag && circ.max<=6*10^mag) circles <- tail(seq(0, 6*10^mag, length.out=5), -1)
		if(circ.max>6*10^mag && circ.max<=8*10^mag) circles <- tail(seq(0, 8*10^mag, length.out=5), -1)
		if(circ.max>8*10^mag) circles <- tail(seq(0, 10^(mag+1), length.out=5), -1)
	}
	
	if(is.numeric(width.leg)) if(width.leg>1) stop("'width.leg' must be a numeric value between 0 and 1")
	if(num.classes>1 && width.leg!=0) lo <- layout(matrix(1:2, 1, 2), widths=c(1, width.leg))
	
	# plot
	par(mar=c(1,1,1,1), las=1)
	plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1] > pin[2]) {
		xlim <- (pin[1]/pin[2]) * xlim
	} else {
		ylim <- (pin[2]/pin[1]) * ylim
	}
	plot.window(xlim, ylim, "", asp=1)
	
	if(!fg) {
		for(c in num.classes:1) {
	    	plot.data <- c(tail(rev(as.vector(e.cum[,c])), n=1), head(rev(as.vector(e.cum[,c])), n=-1))
			
			for (i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.data[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				polygon(xlist, ylist, col=col[c], border=col[c], lwd=0.01)
			}
		}
		if(!is.null(col.border)) {
			e.cum[is.na(e.cum)] <- 0
			plot.max <- apply(e.cum, 1, max)
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
	if(show.total) text(-1, 1, paste("Total:", total), pos=4, cex=cex.axis, col=col.axis)
	if(unit=="kWh/m^2/a") text(1, -1, expression(paste("kWh/m"^2,"/a")), pos=2, cex=cex.axis, col=col.axis)
	else text(1, -1, unit, pos=2, cex=cex.axis, col=col.axis)
	
	if(fg) {
		for(c in num.classes:1) {
	    	plot.data <- c(tail(rev(as.vector(e.cum[,c])), n=1), head(rev(as.vector(e.cum[,c])), n=-1))
			
			for (i in 1:num.sectors) {
				arc.pts <- seq(sector.edges[i] - sector.width/2*sec.space, sector.edges[i] + sector.width/2*sec.space, length.out=trunc(360/num.sectors*sec.space))
				rad <- 0.9 * plot.data[i]/tail(circles, 1)
				xlist <- c(0, rad * cos(arc.pts), 0)
				ylist <- c(0, rad * sin(arc.pts), 0)
				polygon(xlist, ylist, col=col[c], border=col[c], lwd=0.01)
			}
		}
		if(!is.null(col.border)) {
			e.cum[is.na(e.cum)] <- 0
			plot.max <- apply(e.cum, 1, max)
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
		legend("left", legend=names(e.cum), title=title.leg, fill=col[1:num.classes], xjust=0, bty=bty.leg, border=border.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
	}
}
