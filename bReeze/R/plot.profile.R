plot.profile <-
function(x, sector, measured=TRUE, ...) {
###	plotting profile
		
	if(is.null(attr(x, "call")$mast)) stop("Source mast object of ", substitute(x), " could not be found")
	mast <- get(attr(x, "call")$mast)
	v.set <- attr(x, "call")$v.set
	dir.set <- attr(x, "call")$dir.set
	num.sectors <- attr(x, "call")$num.sectors
	subset <- attr(x, "call")$subset
	h.ref <- x$h.ref
	sector.names <- row.names(x$profile)
	
	if(missing(sector)) sector <- NULL
	if(length(sector)>1) stop("Please choose only one 'sector' by name or index")
	if(is.numeric(sector)) {
		if(sector<1 || sector>num.sectors+1) stop("Sector not found")
	} else if(is.character(sector)) {
		sector <- match(sector, sector.names)
		if(is.na(sector)) stop("Sector not found")
	} else {
		if(!is.null(sector)) stop("Sector not found - please choose 'sector' by name or index")
	}
	
	# subset
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) {
		col <- plot.param$col
		if(length(col)==1 && is.null(sector)) col <- rep(col, num.sectors+1)
	} else {
		col <- sample(colors(), num.sectors+1)
		if(num.sectors==4) col <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")
		if(num.sectors==8) col <- c("#377EB8", "#41B6C4", "#4DAF4A", "#9970AB", "#984EA3", "#F781BF", "#FF7F00", "#A6761D", "#E41A1C")
		if(num.sectors==12) col <- c("#08519C", "#3182BD", "#74C476", "#006D2C", "#31A354", "#9E9AC8", "#54278F", "#756BB1", "#FED976", "#FD8D3C", "#FEB24C", "#6BAED6", "#E41A1C")
		if(num.sectors==16) col <- c("#08519C", "#3182BD", "#41B6C4", "#74C476", "#006D2C", "#31A354", "#9970AB", "#9E9AC8", "#54278F", "#756BB1", "#F781BF", "#FED976", "#FD8D3C", "#FEB24C", "#A6761D", "#6BAED6", "#E41A1C")
		if(!is.null(sector)) col <- col[sector]
	}
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "black"
	if(any(names(plot.param)=="col.leg")) col.leg <- plot.param$col.leg
	else col.leg <- "black"
	if(any(names(plot.param)=="col.ticks")) col.ticks <- plot.param$col.ticks
	else col.ticks <- "black"
	if(any(names(plot.param)=="col.box")) col.box <- plot.param$col.box
	else col.box <- "black"
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex
	if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
	else cex.leg <- cex-0.2
	if(any(names(plot.param)=="lty")) {
		lty <- plot.param$lty
		if(length(lty)==1 && is.null(sector)) lty <- rep(lty, num.sectors+1)
	} else {
		lty <- c(rep(5, num.sectors), 1)
		if(num.sectors==4) lty <- c(5, 5, 5, 5, 1)
		if(num.sectors==8) lty <- c(5, 3, 5, 3, 5, 3, 5, 3, 1)
		if(num.sectors==12) lty <- c(5, 4, 3, 5, 4, 3, 5, 4, 3, 5, 4, 3, 1)
		if(num.sectors==16) lty <- c(5, 4, 2, 3, 5, 4, 2, 3, 5, 4, 2, 3, 5, 4, 2, 3, 1)
		if(!is.null(sector)) lty <- lty[sector]
	}	
	if(any(names(plot.param)=="lwd")) {
		lwd <- plot.param$lwd
		if(length(lwd)==1 && is.null(sector)) lwd <- rep(lwd, num.sectors+1)
	} else {
		lwd <- c(rep(1.2, num.sectors), 2)
		if(!is.null(sector)) lwd <- lwd[sector]
	}
	if(any(names(plot.param)=="pch")) pch <- plot.param$pch
	else pch <- 0
	if(any(names(plot.param)=="xlim")) xlim <- plot.param$xlim
	else xlim <- c(0, 1.75*ceiling(max(frequency(mast, v.set=v.set[1], dir.set=dir.set, num.sectors=num.sectors, bins=NULL, subset=subset, print=FALSE)$wind.speed, na.rm=TRUE)))
	if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
	else ylim <- c(0,200)
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
	else pos.leg <- "topright"
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- "Wind speed [m/s]"
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else ylab <- "Height [m]"
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(4,4,1,1)
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(2.2,0.7,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	if(any(names(plot.param)=="legend")) legend <- plot.param$legend
	else legend <- TRUE
	
	# calculate and plot
	h.range <- c(1:max(ylim))
	v.range <- seq(xlim[1],xlim[2],0.1)
	v.mean <- data.frame(matrix(NA, ncol=1, nrow=num.sectors+1))
	h <- NULL
	if(measured) {		
		for(i in 1:length(v.set)) {
			if(!is.null(mast$sets[[v.set[i]]]$data$v.avg[start:end])) {
				v.mean <- data.frame(v.mean, cbind(frequency(mast, v.set=v.set[i], dir.set=dir.set, num.sectors=num.sectors, bins=NULL, subset=subset, print=FALSE)$wind.speed))
				h <- append(h, mast$sets[[v.set[i]]]$height)
				names(v.mean)[i] <- names(mast$sets)[v.set[i]]
			}
		}
	}
	v.mean[1] <- NULL
	
	par(mar=mar, mgp=mgp, las=las, bty="n")
	if(is.null(sector)) { # all sectors
		v.over.h <- x$profile$v.ref[1] * exp(x$profile$alpha[1] * log(h.range / h.ref))
		h.over.v <- spline(x=v.over.h, y=h.range, method="natural", xout=v.range)
		h.over.v[[2]][h.over.v[[2]]<0] <- 0	
		plot(h.over.v, type="l", xlim=xlim, ylim=ylim, axes=FALSE, lty=lty[1], lwd=lwd[1], col=col[1], xlab=xlab, ylab=ylab, cex.lab=cex.lab, col.lab=col.lab)
		box(bty=bty, col=col.box)
		axis(1, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		if(measured) for(j in 1:length(v.mean)) points(x=v.mean[1,j], y=h[j], col=col[1], pch=pch, cex=cex-0.2)
		
		for(i in 2:num.sectors) {
			v.over.h <- x$profile$v.ref[i] * exp(x$profile$alpha[i] * log(h.range / h.ref))
			h.over.v <- spline(x=v.over.h, y=h.range, method="natural", xout=v.range)
			h.over.v[[2]][h.over.v[[2]]<0] <- 0
			lines(h.over.v, lty=lty[i], lwd=lwd[i], col=col[i])
			if(measured) for(j in 1:length(v.mean)) points(x=v.mean[i,j], y=h[j], col=col[i], pch=pch, cex=cex-0.2)
		}
		
		v.over.h <- x$profile$v.ref[num.sectors+1] * exp(x$profile$alpha[num.sectors+1] * log(h.range / h.ref))
		h.over.v <- spline(x=v.over.h, y=h.range, method="natural", xout=v.range)
		h.over.v[[2]][h.over.v[[2]]<0] <- 0
		lines(h.over.v, lty=lty[num.sectors+1], lwd=lwd[num.sectors+1], col=col[num.sectors+1])
		if(measured) for(j in 1:length(v.mean)) points(x=v.mean[num.sectors+1,j], y=h[j], col=col[num.sectors+1], pch=pch, cex=cex-0.2)
		
		if(legend) {
			if(measured) legend(pos.leg, legend=c(sector.names, "measured"), col=c(col,"darkgrey"), lty=c(lty,NA), lwd=c(lwd,NA), pch=c(rep(NA,num.sectors+1),pch), bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
			else legend(pos.leg, legend=sector.names, col=col, lty=lty, lwd=lwd, bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
		}
	} else { # one sector
		v.over.h <- x$profile$v.ref[sector] * exp(x$profile$alpha[sector] * log(h.range / h.ref))
		h.over.v <- spline(x=v.over.h, y=h.range, method="natural", xout=v.range)
		h.over.v[[2]][h.over.v[[2]]<0] <- 0
		plot(h.over.v, type="l", xlim=xlim, ylim=ylim, axes=FALSE, lty=lty, lwd=lwd, col=col, xlab=xlab, ylab=ylab, cex.lab=cex.lab, col.lab=col.lab)
		box(bty=bty, col=col.box)
		axis(1, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		if(measured) for(j in 1:length(v.mean)) points(x=v.mean[sector,j], y=h[j], col=col, pch=pch, cex=cex-0.2)
		
		if(legend) if(measured) legend(pos.leg, legend="measured", col=col, pch=pch, bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, text.col=col.leg)
	}
}
