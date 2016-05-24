plot.month.stats <- 
function(x, set, ...) {
### plotting monthly data
		
	if(is.data.frame(x)) num.sets <- 1
	else num.sets <- length(x)
	if(missing(set)) set <- 1:num.sets
	n.set <- length(set)
	if(!is.numeric(set)) set <- match(set, names(x))
	if(any(is.na(set))) stop("'set' not found")
	if(any(set<1) || any(set>num.sets)) stop("'set' not found")
	unit <- attr(x, "unit")
	years <- length(x[[1]])-2
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) col <- plot.param$col
	else {
		if(n.set<=9) {
			if(requireNamespace("RColorBrewer", quietly=TRUE)) {
				col <- col1 <- RColorBrewer::brewer.pal(3, "Paired")
				if(years>3) col <- col1 <- RColorBrewer::brewer.pal(years, "Paired")
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
	else cex.leg <- cex
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
	else pos.leg <- NULL
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- "Months"
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else ylab <- paste("Wind speed [", unit, "]", sep="")
	if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
	else ylim <- NULL
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(4,5,1,1)
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(2.5,1,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	if(any(names(plot.param)=="col.box")) col.box <- plot.param$col.box
	else col.box <- "black"
	if(any(names(plot.param)=="plot.names")) plot.names <- plot.param$plot.names
	else plot.names <- TRUE
	if(any(names(plot.param)=="legend")) legend <- plot.param$legend
	else legend <- TRUE
	if(any(names(plot.param)=="border")) border <- plot.param$border
	else border <- NA
	
	# plot
	if(n.set==1 || num.sets==1) {
		if(is.null(pos.leg)) pos.leg <- "top"
		par(mar=mar, mgp=mgp, las=las)
		if(is.null(ylim)) ylim <- c(-0.1, ceiling(max(x[[set]][1:12,1:(length(x[[set]])-2)], na.rm=TRUE))+0.3)
		barplot(t(as.matrix(x[[set]][1:12,1:(length(x[[set]])-2)])), beside=TRUE, xaxt="n", yaxt="n", col=col[1:years], border=border, ylim=ylim, xpd=FALSE)
		box(bty=bty, col=col.box)
		axis(2, line=mgp[3], col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		bxp <- barplot(t(as.matrix(x[[set]][1:12,1:(length(x[[set]])-2)])), beside=TRUE, plot=FALSE)
		at <- apply(bxp, 2, mean)
		mtext(toupper(row.names(x[[set]])[1:12]), side=1, line=mgp[2]-0.6, at=at, cex=cex.axis-0.1, col=col.axis)
		mtext(xlab, side=1, line=mgp[1]-0.5, at=mean(at), cex=cex.lab+0.1, col=col.lab, las=1)
		mtext(ylab, side=2, line=mgp[1], las=0, cex=cex.lab+0.1, col=col.lab)
		if(plot.names) mtext(names(x)[set], side=2, line=mgp[1]+1.2, las=0, cex=cex.lab+0.1, col=col.lab)
		if(legend) legend(pos.leg, legend=names(x[[1]])[1:years], fill=col[1:years], border=border, ncol=years, bty=bty.leg, cex=cex.leg-0.1, x.intersp=x.intersp, text.col=col.leg)
	} else {
		if(is.null(pos.leg)) pos.leg <- "center"
		lo <- layout(matrix(c(n.set+2, 1:(n.set+1)), n.set+2, 1), heights=c(1, rep(5, n.set), 1))
		par(mar=c(1,5.5,0,1), mgp=mgp, las=las)
		dat.max <- ceiling(max(unlist(x), na.rm=TRUE))
		for(i in 1:n.set) {
			if(is.null(ylim)) ylim <- c(-0.1, dat.max+0.3)
			barplot(t(as.matrix(x[[i]][1:12,1:years])), beside=TRUE, xaxt="n", yaxt="n", col=col[1:years], border=border, ylim=ylim, xpd=FALSE)
			box(bty=bty, col=col.box)
			axis(2, line=mgp[3], col=col.ticks, col.axis=col.axis, cex.axis=cex.axis+0.2)
			if(plot.names) mtext(names(x)[set[i]], side=2, line=mgp[1]+1.2, las=0, cex=cex.lab, col=col.lab)
			mtext(ylab, side=2, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
		}
		bxp <- barplot(t(as.matrix(x[[1]][1:12,1:(length(x[[1]])-2)])), beside=TRUE, plot=FALSE)
		at <- apply(bxp, 2, mean)
		mtext(toupper(row.names(x[[1]])[1:12]), side=1, line=mgp[2]-0.5, at=at, cex=cex.axis-0.2, col=col.axis)
		mtext(xlab, side=1, line=mgp[1]-0.4, at=mean(at), cex=cex.lab, col=col.lab, las=1)
		plot(0, type="n", axes=FALSE, xlab="", ylab="")
		par(mar=c(0,5.5,0,1))
		plot(0, type="n", axes=FALSE, xlab="", ylab="")
		if(legend) legend(pos.leg, legend=names(x[[1]])[1:years], fill=col[1:years], border=border, ncol=years, bty=bty.leg, cex=cex.leg+0.2, x.intersp=x.intersp, text.col=col.leg)
	}
}
