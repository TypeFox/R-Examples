turb.iec.plot <-
function(mast, set, subset, ...) {
### plotting turbulence intesity and site classification after IEC from mast object
		
	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	num.sets <- length(mast$sets)
	if(!is.numeric(set)) set <- match(set, names(mast$sets))
	if(is.na(set)) stop("'set' not found")
	if(set<0 || set>num.sets) stop("'set' not found")
	if(is.null(mast$sets[[set]]$data$turb.int)) stop("'set' does not contain turbulence intensity data")
	unit <- attr(mast$sets[[set]]$data$v.avg, "unit")
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	vmax <- ceiling(max(mast$sets[[set]]$data$v.avg[start:end], na.rm=TRUE))
	site.turb <- c()
	for(i in 0:(vmax-1)) {
		site.turb <- append(site.turb, mean(mast$sets[[set]]$data$turb.int[mast$sets[[set]]$data$v.avg[start:end]>=i & mast$sets[[set]]$data$v.avg[start:end]<i+1], na.rm=TRUE))
	}
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else col <- "#E41A1C"
	if(any(names(plot.param)=="line")) line <- plot.param$line
	else line <- "black"
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
	if(any(names(plot.param)=="border")) border <- plot.param$border
	else border <- col
	if(any(names(plot.param)=="space")) {
		if(plot.param$space<1 && plot.param$space>0) space <- plot.param$space
		else space <- 0.2
	} else space <- 0.2
	if(any(names(plot.param)=="lty")) lty <- plot.param$lty
	else lty <- c(3, 2, 1)
	if(any(names(plot.param)=="lwd")) lwd <- plot.param$lwd
	else lwd <- 1.2
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex
	if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
	else cex.leg <- cex-0.2
	if(any(names(plot.param)=="xlim")) xlim <- plot.param$xlim
	else xlim <- c(0, vmax)
	if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
	else ylim <- c(0, 0.6)
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
	else pos.leg <- "topright"
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- paste("Wind speed [", unit, "]", sep="")
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else ylab <- "Turbulence intensity [-]"
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(4.5,4.5,1,1)
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(2.2,0.7,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	if(any(names(plot.param)=="legend")) legend <- plot.param$legend
	else legend <- TRUE
	if(any(names(plot.param)=="leg.text")) leg.text <- plot.param$leg.text
	else leg.text <- c("Class A (0.16)", "Class B (0.14)", "Class C (0.12)", "Site")
	
	if(length(line)==1) line <- rep(line, 3)
	if(length(lty)==1) lty <- rep(lty, 3)
	if(length(lwd)==1) lwd <- rep(lwd, 3)
	
	v <- seq(0, xlim[2], 1)
	sigma1 <- 0.16*(0.75*v+5.6)/v
	sigma2 <- 0.14*(0.75*v+5.6)/v
	sigma3 <- 0.12*(0.75*v+5.6)/v
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	par(mar=mar, mgp=mgp, las=las, bty="n")
	
	# plot
	plot(v, sigma1, type="l",  xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, axes=FALSE, lty=lty[3], lwd=lwd[3], col=line[3], cex.lab=cex.lab, col.lab=col.lab)
	box(bty=bty, col=col.box)
	axis(1, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
	axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
	lines(v, sigma2, lty=lty[2], lwd=lwd[2], col=line[2])
	lines(v, sigma3, lty=lty[1], lwd=lwd[1], col=line[1])
	for(i in 5:vmax) {
		polygon(c(i-space/2, i-space/2, i-1+space/2, i-1+space/2), c(0, site.turb[i], site.turb[i], 0), col=col, border=border)
	}
	if(legend) legend(pos.leg, legend=leg.text, col=c(line, border), lty=c(lty, NA), lwd=c(lwd, NA), pch=c(NA, NA, NA, 22), pt.bg=c(NA, NA, NA, col), bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
}
