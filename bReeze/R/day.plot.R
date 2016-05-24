day.plot <-
function(mast, set, dir.set=set, signal, num.sectors=NULL, subset, ...) {
### plotting diurnal wind speed data
	
	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	num.sets <- length(mast$sets)
	if(missing(set)) set <- "all"
	if(missing(signal)) stop("No signal to plot")
	if(length(signal)>1) stop("Please choose only one signal")
	if(!is.null(num.sectors)) {
		if(is.null(dir.set)) stop("Sectoral plot requires dir.set")
		if(!is.numeric(num.sectors)) stop("'num.sectors' must be numeric or NULL")
		if(num.sectors<=1) stop("'num.sectors' must be greater 1")
	}
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
		
	h.unit <- attr(mast$sets[[1]]$height, "unit")
	unit <- NULL
	for(i in 1:num.sets) {
		if(any(names(mast$sets[[i]]$data)==signal)) unit <- attr(mast$sets[[i]]$data[,signal], "unit")
		break
	}
	
	# prepare plot
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) {
		col <- plot.param$col
		if(!is.null(num.sectors) && length(col)==1) col <- rep(col, num.sectors+1)
	} else {
		if(is.null(num.sectors)) {
			if(num.sets<=9) {
				if(suppressWarnings(require(RColorBrewer, quietly=TRUE))) {
					col <- col1 <- RColorBrewer::brewer.pal(3, "Set1")
					if(num.sets>3) col <- col1 <- RColorBrewer::brewer.pal(num.sets, "Set1")
					col[1] <- col1[2]
					col[2] <- col1[1]
				} else col <- c("blue", "green", "red", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
			} else col <- c("blue", "green", "red", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
		} else {
			col <- sample(colors(), num.sectors+1)
			if(num.sectors==4) col <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")
			if(num.sectors==8) col <- c("#377EB8", "#41B6C4", "#4DAF4A", "#9970AB", "#984EA3", "#F781BF", "#FF7F00", "#A6761D", "#E41A1C")
			if(num.sectors==12) col <- c("#08519C", "#3182BD", "#74C476", "#006D2C", "#31A354", "#9E9AC8", "#54278F", "#756BB1", "#FED976", "#FD8D3C", "#FEB24C", "#6BAED6", "#E41A1C")
			if(num.sectors==16) col <- c("#08519C", "#3182BD", "#41B6C4", "#74C476", "#006D2C", "#31A354", "#9970AB", "#9E9AC8", "#54278F", "#756BB1", "#F781BF", "#FED976", "#FD8D3C", "#FEB24C", "#A6761D", "#6BAED6", "#E41A1C")
		}
	}
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "black"
	if(any(names(plot.param)=="col.ticks")) col.ticks <- plot.param$col.ticks
	else col.ticks <- "black"
	if(any(names(plot.param)=="col.box")) col.box <- plot.param$col.box
	else col.box <- "black"
	if(any(names(plot.param)=="col.leg")) col.leg <- plot.param$col.leg
	else col.leg <- "black"
	if(any(names(plot.param)=="lty")) {
		lty <- plot.param$lty
		if(!is.null(num.sectors) && length(lty)==1) lty <- rep(lty, num.sectors+1)
	} else {
		if(is.null(num.sectors)) lty <- rep(1, num.sets)
		else {
			lty <- c(rep(5, num.sectors), 1)
			if(num.sectors==4) lty <- c(5, 5, 5, 5, 1)
			if(num.sectors==8) lty <- c(5, 3, 5, 3, 5, 3, 5, 3, 1)
			if(num.sectors==12) lty <- c(5, 4, 3, 5, 4, 3, 5, 4, 3, 5, 4, 3, 1)
			if(num.sectors==16) lty <- c(5, 4, 2, 3, 5, 4, 2, 3, 5, 4, 2, 3, 5, 4, 2, 3, 1)
		}
	}
	if(any(names(plot.param)=="lwd")) {
		lwd <- plot.param$lwd
		if(!is.null(num.sectors) && length(lwd)==1) lwd <- rep(lwd, num.sectors+1)
	} else {
		if(is.null(num.sectors)) lwd <- rep(1, num.sets)
		else lwd <- c(rep(1, num.sectors), 2)
	}
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex
	if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
	else cex.leg <- cex-0.2
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
	else pos.leg <- "topright"
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- "Hour of day"
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else {
		ylab <- signal
		if(signal=="v.avg" || signal=="v.max" || signal=="v.min") ylab <- paste0("Wind speed [", unit, "]")
		if(signal=="dir.avg") ylab <- paste0("Wind direction [", unit, "]")
		if(signal=="turb.int") ylab <- paste0("Turbulence intensity [", unit, "]")
	}
	if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
	else ylim <- NULL
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(4.5,5,1.5,1)
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(2,0.7,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	par(mar=mar, mgp=mgp, las=las, bty="n")
	
	# calculate and plot	
	if(set!="all") { # one set
		if(!is.numeric(set)) if(set!="all") set <- match(set, names(mast$sets))
		if(is.na(set)) stop("'set' not found")
		if(set<0 || set>num.sets) stop("'set' not found")
		if(!any(names(mast$sets[[set]]$data)==signal)) stop("'set' does not contain the choosen signal")
		dat <- mast$sets[[set]]$data[,which(names(mast$sets[[set]]$data)==signal)][start:end]
		if(!is.numeric(dir.set)) dir.set <- match(dir.set, names(mast$sets))
		if(is.na(dir.set)) stop("'dir.set' not found")
		if(dir.set<0 || dir.set>num.sets) stop("'dir.set' not found")
		if(!any(names(mast$sets[[dir.set]]$data)=="dir.avg")) stop("'dir.set' does not contain wind direction data")
		if(!is.null(num.sectors)) {
			sector.width <- 360/num.sectors
			sectors <- seq(0, 360-sector.width, by=sector.width)
			sector.edges <- c(sectors-sector.width/2, tail(sectors, n=1)+sector.width/2)%%360
		}
		if(length(col)==1) col <- rep(col, set)
		if(length(lty)==1) lty <- rep(lty, set)
		if(length(lwd)==1) lwd <- rep(lwd, set)
		
		# mean diurnal
		diurnal <- NULL
		for(i in 0:23) {
			hour.idx <- mast$timestamp[start:end]$hour==i
			hour.values <- dat[hour.idx]
			hour.values <- hour.values[!is.na(hour.values)]
			if(length(hour.values)>0) diurnal <- append(diurnal, mean(hour.values))
		}
		diurnal <- list(diurnal)
		
		if(!is.null(num.sectors)) { # sectoral
			for(sec in 1:num.sectors) {
				low <- sector.edges[sec]
				high <- sector.edges[sec+1]
				if(low<high) idx.dir <- mast$sets[[dir.set]]$data$dir.avg[start:end]>=low & mast$sets[[dir.set]]$data$dir.avg[start:end]<high
				else idx.dir <- mast$sets[[dir.set]]$data$dir.avg[start:end]>=low | mast$sets[[dir.set]]$data$dir.avg[start:end]<high
				diurnal.s <- NULL
				for(i in 0:23) {
					hour.idx <- mast$timestamp[start:end]$hour==i
					hour.values <- dat[hour.idx & idx.dir]
					hour.values <- hour.values[!is.na(hour.values)]
					if(length(hour.values)>0) diurnal.s <- append(diurnal.s, mean(hour.values))
				}
				diurnal[[length(diurnal)+1]] <- diurnal.s
			}
		}
		if(is.null(ylim)) {
			ylim <- c(0.8*min(diurnal[[1]]), 1.2*max(diurnal[[1]]))
			if(!is.null(num.sectors)) for(d in 2:length(diurnal)) ylim <- c(min(ylim[1], 0.8*min(diurnal[[d]])), max(ylim[2], 1.2*max(diurnal[[d]])))
		}
		if(length(diurnal)==1) plot(0:23, diurnal[[1]], type="l", xaxt="n", yaxt="n", xlim=c(0, 24), ylim=ylim, xlab=xlab, ylab=ylab, col=col[set], lty=lty[set], lwd=lwd[set], cex.lab=cex.lab, col.lab=col.lab)
		else plot(0:23, diurnal[[1]], type="n", xaxt="n", yaxt="n", xlim=c(0, 24), ylim=ylim, xlab=xlab, ylab=ylab, cex.lab=cex.lab, col.lab=col.lab)
		box(bty=bty, col=col.box)
		axis(1, at=c(0,6,12,18,24), col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)

		if(length(diurnal)>1) {
			for(d in 2:length(diurnal)) lines(0:23, diurnal[[d]], col=col[d-1], lty=lty[d-1], lwd=lwd[d-1])
			lines(0:23, diurnal[[1]], col=col[length(diurnal)], lty=lty[length(diurnal)], lwd=lwd[length(diurnal)])
		}
		
		#mtext(xlab, 1, 2, cex=cex.lab, col=col.lab)
		if(!is.null(pos.leg)) {
			if(!is.null(num.sectors)) {
				sec <- c(paste0("s", 1:num.sectors),"all")
				if(num.sectors==4) sec <- c("n","e","s","w","all")
				if(num.sectors==8) sec <- c("n","ne","e","se","s","sw","w","nw","all")
				if(num.sectors==12) sec <- c("n","nne","ene","e","ese","sse","s","ssw","wsw","w","wnw","nnw","all")
				if(num.sectors==16) sec <- c("n","nne","ne","ene","e","ese","se","sse","s","ssw","sw","wsw","w","wnw","nw","nnw","all")
				legend(pos.leg, legend=sec, col=col, lty=lty, lwd=lwd, bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
			} else legend(pos.leg, legend=paste0(names(mast$sets)[set], " (", mast$sets[[set]]$height, h.unit, ")"), col=col[set], lty=lty[set], lwd=lwd[set], bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
		}
	} else { # all sets
		if(!is.null(num.sectors)) stop("Sectoral plot not available for multiple sets")
		set.index <- NULL
		for(s in 1:num.sets) if(any(names(mast$sets[[s]]$data)==signal)) set.index <- append(set.index, s)
		if(is.null(set.index)) stop("Signal not found in any set")
		if(any(names(plot.param)=="col")) {
			n.set <- length(set.index)
			if(length(col)==1) col <- rep(col, n.set)
			if(n.set!=length(col)) stop(n.set, " colours needed")
			set.all <- 1:set.index[n.set]
			col.all <- rep(NA, set.index[n.set])
			col.all[set.index] <- col
			col <- col.all
		}
		if(any(names(plot.param)=="lty")) {
			n.set <- length(set.index)
			if(length(lty)==1) lty <- rep(lty, n.set)
			if(n.set!=length(lty)) stop(n.set, " line types needed")
			set.all <- 1:set.index[n.set]
			lty.all <- rep(NA, set.index[n.set])
			lty.all[set.index] <- lty
			lty <- lty.all
		}
		if(any(names(plot.param)=="lwd")) {
			n.set <- length(set.index)
			if(length(lwd)==1) lwd <- rep(lwd, n.set)
			if(n.set!=length(lwd)) stop(n.set, " line widths needed")
			set.all <- 1:set.index[n.set]
			lwd.all <- rep(NA, set.index[n.set])
			lwd.all[set.index] <- lwd
			lwd <- lwd.all
		}
		
		diurnal <- NULL
		for(i in 0:23) {
			hour.idx = mast$timestamp[start:end]$hour==i
			hour.values <- mast$sets[[set.index[1]]]$data[hour.idx,which(names(mast$sets[[set.index[1]]]$data)==signal)]
			hour.values <- hour.values[!is.na(hour.values)]
			if(length(hour.values)>0) diurnal <- append(diurnal, mean(hour.values))
		}
		if(is.null(ylim)) ylim <- c(0.8*min(diurnal), 1.2*max(diurnal))
		plot(0:23, diurnal, type="l", xaxt="n", yaxt="n", xlim=c(0, 24), ylim=ylim, xlab=xlab, ylab=ylab, col=col[set.index[1]], lty=lty[set.index[1]], lwd=lwd[set.index[1]], cex.lab=cex.lab, col.lab=col.lab)
		box(bty=bty, col=col.box)
		axis(1, at=c(0,6,12,18,24), col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)

		if(length(set.index)>1) {
			for(s in 2:length(set.index)) {
				diurnal <- NULL
				for(i in 0:23) {
					hour.idx <- mast$timestamp[start:end]$hour==i
					hour.values <- mast$sets[[set.index[s]]]$data[hour.idx,which(names(mast$sets[[set.index[s]]]$data)==signal)]
					hour.values <- hour.values[!is.na(hour.values)]
					if(length(hour.values)>0) diurnal <- append(diurnal, mean(hour.values))
				}
				lines(0:23, diurnal, col=col[set.index[s]], lty=lty[set.index[s]], lwd=lwd[set.index[s]])
			}
		}
		
		heights <- c()
		for(s in 1:num.sets) {
			if(any(names(mast$sets[[s]]$data)==signal)	) heights <- append(heights, mast$sets[[s]]$height)
		}
		
		if(!is.null(pos.leg)) legend(pos.leg, legend=paste0(names(mast$sets)[set.index], " (", heights, h.unit, ")"), col=col[set.index], lty=lty[set.index], lwd=lwd[set.index], bty=bty.leg, cex=cex.leg, x.intersp=x.intersp, y.intersp=y.intersp, text.col=col.leg)
	}
}
