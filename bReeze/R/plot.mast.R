plot.mast <-
function(x, set, signal=c("v.avg", "dir.avg", "turb.int"), subset, ...) {
### plotting time series of mast data
		
	num.sets <- length(x$sets)
	timestamp <- x$timestamp
	
	if(missing(set)) set <- 1:num.sets
	if(!is.numeric(set)) set <- match(set, names(x$sets))
	if(any(is.na(set))) stop("'set' not found\n")
	if(any(set<1) || any(set>num.sets)) stop("'set' not found\n")
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	# get units
	n.sig <- length(signal)
	n.set <- length(set)
	h.unit <- attr(x$sets[[1]]$height, "unit")
	units <- rep("", n.sig)
	for(s in 1:n.sig) {
		for(i in 1:n.set) {
			if(any(names(x$sets[[i]]$data)==signal[s])) {
				if(!is.null(attr(x$sets[[i]]$data[,signal[s]], "unit")))	units[s] <- attr(x$sets[[i]]$data[,signal[s]], "unit"); break
			}
		}
	}
	  
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else {
		if(num.sets<=9) {
			if(requireNamespace("RColorBrewer", quietly=TRUE)) {
				col <- col1 <- RColorBrewer::brewer.pal(3, "Set1")
				if(num.sets>3) col <- col1 <- RColorBrewer::brewer.pal(num.sets, "Set1")
				col[1] <- col1[2]
				col[2] <- col1[1]
			} else col <- c("blue", "green", "red", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
		} else col <- c("blue", "green", "red", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
	}
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "black"
	if(any(names(plot.param)=="col.leg")) col.leg <- plot.param$col.leg
	else col.leg <- "black"
	if(any(names(plot.param)=="col.ticks")) col.ticks <- plot.param$col.ticks
	else col.ticks <- "black"
	if(any(names(plot.param)=="cex")) cex <- plot.param$cex
	else cex <- 1
	cex <- cex-0.2
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- cex
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- cex
	if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
	else cex.leg <- cex-0.1
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(1,4.5,0,1)
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(2.5,0.7,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	if(any(names(plot.param)=="col.box")) col.box <- plot.param$col.box
	else col.box <- "black"
	if(any(names(plot.param)=="legend")) legend <- plot.param$legend
	else legend <- TRUE
	if(any(names(plot.param)=="lty")) lty <- plot.param$lty
	else lty <- rep(1, num.sets)
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else {
		ylab <- NULL
		for(i in 1:n.sig) {
			if(signal[i]=="v.avg") ylab <- append(ylab, paste0("Wind speed [", units[i], "]"))
			else if(signal[i]=="v.max") ylab <- append(ylab, paste0("Max wind speed [", units[i], "]"))
			else if(signal[i]=="v.min") ylab <- append(ylab, paste0("Min wind speed [", units[i], "]"))
			else if(signal[i]=="dir.avg") ylab <- append(ylab, paste0("Wind direction [", units[i], "]"))
			else if(signal[i]=="turb.int") ylab <- append(ylab, paste0("Turbulence intensity [", units[i], "]"))
			else ylab <- append(ylab, paste0(signal[i], " [", units[i], "]"))
		}
	}
	for(i in 1:length(ylab)) if(substr(ylab[i], nchar(ylab[i])-2, nchar(ylab[i]))==" []") ylab[i] <- substr(ylab[i], 1, nchar(ylab[i])-3)
	
	# layout
	lo <- layout(matrix(c(n.sig+2, 1:(n.sig+1)), n.sig+2, 1), heights=c(1, rep(4, n.sig), 1))
	if(n.sig==1) lo <- layout(matrix(c(n.sig+2, 1:(n.sig+1)), n.sig+2, 1), heights=c(1, 9, 1))
	n.set <- length(set)
	set.idx <- data.frame(matrix(NA, ncol=n.sig, nrow=n.set))
	names(set.idx) <- signal
	for(i in 1:n.sig) {
		for(j in 1:n.set) {
			if(any(names(x$sets[[set[j]]]$data)==signal[i])) set.idx[j,i] <- j
		}
	}
		
	# plot
	for(i in 1:n.sig) {
		par(mar=mar, mgp=mgp, las=las, bty="n")
		sets <- set.idx[!is.na(set.idx[,which(names(set.idx)==signal[i])]),which(names(set.idx)==signal[i])]
		if(length(sets)>=1) {
			plot(timestamp[start:end], x$sets[[sets[1]]]$data[[which(names(x$sets[[sets[1]]]$data)==signal[i])]][start:end], type="l", col=col[sets[1]], ylab=ylab[i], axes=FALSE, col.lab=col.lab, cex.lab=cex.lab, lty=lty[sets[1]])
			box(bty=bty, col=col.box)
			axis(2, line=mgp[3], col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
			if(i<n.sig) axis.POSIXct(1, at=seq(min(timestamp[start:end]), max(timestamp[start:end]), length.out=6), format="%Y-%m-%d %H:%M:%S", labels=FALSE, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
			else axis.POSIXct(1, at=seq(min(timestamp[start:end]), max(timestamp[start:end]), length.out=6), format="%Y-%m-%d %H:%M:%S", col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		
		
			if(length(sets)>1) {
				for(j in 2:length(sets)) {
					lines(timestamp[start:end], x$sets[[sets[j]]]$data[[which(names(x$sets[[sets[j]]]$data)==signal[i])]][start:end], col=col[sets[j]], lty=lty[sets[j]])
				}
			}
		} else {
			plot(0, type="n", axes=FALSE, xlab="", ylab="")
			text(0, labels=paste(signal[i], "not found!"))
		}
	}
	
	set.idx <- unique(unlist(set.idx)[!is.na(unlist(set.idx))])
	heights <- names <- NULL
	for(i in 1:length(set.idx)) {
		heights <- append(heights, x$sets[[set.idx[i]]]$height)
		names <- append(names, names(x$sets)[set.idx[i]])
	}
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	par(mar=c(0,5,0,1))
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	if(legend) legend("center", legend=paste0(names, " (", heights, h.unit, ")"), col=col[set.idx], lty=lty[set.idx], ncol=length(set.idx), bty=bty.leg, cex=cex.leg, text.col=col.leg, x.intersp=x.intersp)
}
