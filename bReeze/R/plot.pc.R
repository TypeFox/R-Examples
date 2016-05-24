plot.pc <-
function(x, cp=TRUE, ct=TRUE, ...) {
###	plotting power curve object
	
	#pc.list <- FALSE
	#if(class(x)!="pc") {
	#	if(!is.list(x)) stop(substitute(x), " is no power curve object - use createPC to create a power curve or readPC to import from file")
	#	else {
	#		pc.list <- TRUE
	#		for(i in 1:length(x)) {
	#			if(attr(x[[i]], "call")$func!="createPC" && attr(x[[i]], "call")$func!="readPC") stop("At least one object of ", substitute(x), " is not a power curve object")
	#		}
	#	}
	#}
		
	#if(!pc.list) unit <- attr(x, "units")
	#else unit <- attr(x[[1]], "units")
	unit <- attr(x, "units")
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	#if(!pc.list) {
		if(any(names(plot.param)=="col")) col <- plot.param$col
		else col <- c("#3182BD", "#E41A1C", "#E41A1C")
	#} else {
	#	if(any(names(plot.param)=="col")) col <- plot.param$col
	#	else {
	#		if(length(x)<=9) {
	#			if(suppressWarnings(require(RColorBrewer, quietly=TRUE))) {
	#				col <- brewer.pal(3, "Set1")
	#				if(length(x)==2) col <- col[1:2]
	#				if(length(x)>3) col <- brewer.pal(length(x), "Set1")
	#			} else col <- c("blue", "green", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
	#		} else col <- c("blue", "green", "cyan", "magenta", "orange", "brown", "violet", "yellow", "pink", colors())
	#	}
	#}
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
	if(any(names(plot.param)=="lty")) lty <- plot.param$lty
	else lty <- c(1:3)
	if(any(names(plot.param)=="lwd")) lwd <- plot.param$lwd
	else lwd <- c(1,1,1)
	if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
	else x.intersp <- 0.4
	if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
	else y.intersp <- 0.8
	if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
	else bty.leg <- "n"
	if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
	else pos.leg <- "topleft"
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- paste0("Wind speed [", unit[1], "]")
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else ylab <- paste0("Power [", unit[2], "]")
	if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
	else ylim <- NULL
	if(any(names(plot.param)=="xlim")) xlim <- plot.param$xlim
	else xlim <- NULL
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- NULL
	if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
	else mgp <- c(3,1,0)
	if(any(names(plot.param)=="las")) las <- plot.param$las
	else las <- 1
	if(any(names(plot.param)=="bty")) bty <- plot.param$bty
	else bty <- "o"
	if(any(names(plot.param)=="legend")) legend <- plot.param$legend
	else legend <- TRUE
	if(any(names(plot.param)=="leg.text")) leg.text <- plot.param$leg.text
	else leg.text <- NULL
	
	# plot
	#if(!pc.list) {	
		if((cp && !is.null(x$cp)) || (ct && !is.null(x$ct))) {
			if(is.null(mar)) par(mar=c(5,5,1,5), mgp=mgp, las=las, bty="n") else par(mar=mar, mgp=mgp, las=las, bty="n")
		} else {
			if(is.null(mar)) par(mar=c(5,5,1,1), mgp=mgp, las=las, bty="n") else par(mar=mar, mgp=mgp, las=las, bty="n")
		}
		
		plot(x$v[!is.na(x$P)], x$P[!is.na(x$P)], type="l", xaxt="n", yaxt="n", xlab=xlab, ylab=ylab[1], col=col[1], lty=lty[1], lwd=lwd[1], cex=cex, xlim=xlim, ylim=ylim, col.axis=col.axis, col.lab=col.lab, bty="n")
		box(bty=bty, col=col.box)
		axis(1, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		
		if(cp && !is.null(x$cp)) {
			par(new=TRUE, mgp=mgp, las=las)
			plot(x$v[!is.na(x$cp)], x$cp[!is.na(x$cp)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), , col=col[2], lty=lty[2], lwd=lwd[2]) 
			axis(4, col=col.ticks, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, col.axis=col.axis, col.lab=col.lab)  
		}
		if(cp && !is.null(x$cp) && ct && !is.null(x$ct)) {
			par(new=TRUE, mgp=mgp, las=las)
			plot(x$v[!is.na(x$ct)], x$ct[!is.na(x$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), , col=col[3], lty=lty[3], lwd=lwd[3], cex=cex, cex.axis=cex, cex.lab=cex) 
		}
		if((!cp || is.null(x$cp)) && ct && !is.null(x$ct)) {
			par(new=TRUE, mgp=mgp, las=las)
			plot(x$v[!is.na(x$ct)], x$ct[!is.na(x$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), , col=col[2], lty=lty[2], lwd=lwd[2], cex=cex, cex.axis=cex, cex.lab=cex) 
			axis(4, col=col.ticks, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, col.axis=col.axis, col.lab=col.lab)  
		}
		
		if(cp && !is.null(x$cp) && ct && !is.null(x$ct)) {
			if(length(ylab)==1) txt <- "Coefficients [-]" else txt <- ylab[2]
			if(is.null(leg.text)) leg.text <- c("Power", expression(c[p]), expression(c[t]))
			mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
			if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=col, lty=lty[1:3], lwd=lwd[1:3], x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
		} else if(cp && !is.null(x$cp) && (!ct || is.null(x$ct))) {
			if(length(ylab)==1) txt <- "Power coefficient [-]" else txt <- ylab[2]
			if(is.null(leg.text)) leg.text <- c("Power", expression(c[p]))
			mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
			if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=col[c(1,2)], lty=lty[c(1,2)], lwd=lwd[c(1,2)], x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
		} else if(ct && !is.null(x$ct) && (!cp || is.null(x$cp))) {
			if(length(ylab)==1) txt <- "Thrust coefficient [-]" else txt <- ylab[2]
			if(is.null(leg.text)) leg.text <- c("Power", expression(c[t]))
			mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
			if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=col[c(1,2)], lty=lty[c(1,2)], lwd=lwd[c(1,2)], x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
		}
	#} else { # list of pc
	#	cpt <- FALSE
	#	if((cp && !is.null(pc[[1]]$cp)) || (ct && !is.null(pc[[1]]$ct))) cpt <- TRUE
	#	v.min <- min(pc[[1]]$v, na.rm=TRUE)
	#	v.max <- max(pc[[1]]$v, na.rm=TRUE)
	#	P.max <- max(pc[[1]]$P, na.rm=TRUE)
	#	
	#	for(i in 2:length(pc)) {
	#		if((cp && !is.null(pc[[i]]$cp)) || (ct && !is.null(pc[[i]]$ct))) cpt <- TRUE
	#		if(min(pc[[i]]$v, na.rm=TRUE)<v.min) v.min <- min(pc[[i]]$v, na.rm=TRUE)
	#		if(max(pc[[i]]$v, na.rm=TRUE)>v.max) v.max <- max(pc[[i]]$v, na.rm=TRUE)
	#		if(max(pc[[i]]$P, na.rm=TRUE)>P.max) P.max <- max(pc[[i]]$P, na.rm=TRUE)
	#	}	
	#	
	#	if(cpt) {
	#		if(is.null(mar)) par(mar=c(5,5,1,5), mgp=mgp, las=las, bty="n") else par(mar=mar, mgp=mgp, las=las, bty="n")
	#	} else {
	#		if(is.null(mar)) par(mar=c(5,5,1,1), mgp=mgp, las=las, bty="n") else par(mar=mar, mgp=mgp, las=las, bty="n")
	#	}
	#	if(is.null(ylim)) ylim <- c(0, round(P.max))
	#	if(is.null(xlim)) xlim <- c(round(v.min), round(v.max))
	#	
	#	# power curve
	#	plot(pc[[1]]$v[!is.na(pc[[1]]$P)], pc[[1]]$P[!is.na(pc[[1]]$P)], type="l", xaxt="n", yaxt="n", xlab=xlab, ylab=ylab[1], col=col[1], lty=lty[1], lwd=lwd[1], cex=cex, xlim=xlim, ylim=ylim, col.axis=col.axis, col.lab=col.lab, bty="n")
	#	for(i in 2:length(pc)) lines(pc[[i]]$v[!is.na(pc[[i]]$P)], pc[[i]]$P[!is.na(pc[[i]]$P)], col=col[i], lty=lty[1], lwd=lwd[1])
	#	box(bty=bty, col=col.box)
	#	axis(1, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
	#	axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
	#	
	#	# coefficients
	#	if(cp && !is.null(pc[[1]]$cp)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[1]]$v[!is.na(pc[[1]]$cp)], pc[[1]]$cp[!is.na(pc[[1]]$cp)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[1], lty=lty[2], lwd=lwd[2])
	#	}
	#	for(i in 2:length(pc)) if(cp && !is.null(pc[[i]]$cp)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[i]]$v[!is.na(pc[[i]]$cp)], pc[[i]]$cp[!is.na(pc[[i]]$cp)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[i], lty=lty[2], lwd=lwd[2])
	#	}
	#	
	#	if(cp && !is.null(pc[[1]]$cp) && ct && !is.null(pc[[1]]$ct)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[1]]$v[!is.na(pc[[1]]$ct)], pc[[1]]$ct[!is.na(pc[[1]]$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[1], lty=lty[3], lwd=lwd[3], cex=cex, cex.axis=cex, cex.lab=cex)
	#	}
	#	for(i in 2:length(pc)) if(cp && !is.null(pc[[i]]$cp) && ct && !is.null(pc[[i]]$ct)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[i]]$v[!is.na(pc[[i]]$ct)], pc[[i]]$ct[!is.na(pc[[i]]$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[i], lty=lty[3], lwd=lwd[3], cex=cex, cex.axis=cex, cex.lab=cex)
	#	}
	#	
	#	if((!cp || is.null(pc[[1]]$cp)) && ct && !is.null(pc[[1]]$ct)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[1]]$v[!is.na(pc[[1]]$ct)], pc[[1]]$ct[!is.na(pc[[1]]$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[1], lty=lty[2], lwd=lwd[2], cex=cex, cex.axis=cex, cex.lab=cex)
	#	}
	#	for(i in 2:length(pc)) if((!cp || is.null(pc[[i]]$cp)) && ct && !is.null(pc[[i]]$ct)) {
	#		par(new=TRUE, mgp=mgp, las=las)
	#		plot(pc[[i]]$v[!is.na(pc[[i]]$ct)], pc[[i]]$ct[!is.na(pc[[i]]$ct)], type="l", axes=FALSE, ylab="", xlab="", xlim=xlim, ylim=c(0,1), col=col[i], lty=lty[2], lwd=lwd[2], cex=cex, cex.axis=cex, cex.lab=cex)
	#	}
	#	
	#	if(cpt) axis(4, col=col.ticks, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, col.axis=col.axis, col.lab=col.lab)	
	#
	#	# label and legend
	#	cp.p <- ct.p <- FALSE
	#	if(cp && !is.null(pc[[1]]$cp)) cp.p <- TRUE
	#	if(ct && !is.null(pc[[1]]$ct)) ct.p <- TRUE
	#	if(is.null(names(pc))) pc.names <- paste("Power curve", 1:length(pc))
	#	else pc.names <- names(pc)
	#	
	#	for(i in 2:length(pc)) {
	#		if(cp && !is.null(pc[[i]]$cp)) cp.p <- TRUE
	#		if(ct && !is.null(pc[[i]]$ct)) ct.p <- TRUE
	#	}
	#	if(cp.p && ct.p) {
	#		if(length(ylab)==1) txt <- "Coefficients [-]" else txt <- ylab[2]
	#		if(is.null(leg.text)) leg.text <- c(pc.names, NA, "Power", expression(c[p]), expression(c[t]))
	#		mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
	#		if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=c(col, NA, rep("darkgray",3)), lty=c(rep(lty[1], length(pc)), NA, lty[1:3]), lwd=c(rep(lwd[1], length(pc)), NA, lwd[1:3]), x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
	#	} else if(cp.p && !ct.p) {
	#		if(length(ylab)==1) txt <- "Power coefficient [-]" else txt <- ylab[2]
	#		if(is.null(leg.text)) leg.text <- c(pc.names, NA, "Power", expression(c[p]))
	#		mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
	#		if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=c(col, NA, rep("darkgray",2)), lty=c(rep(lty[1], length(pc)), NA, lty[c(1,2)]), lwd=c(rep(lwd[1], length(pc)), NA, lwd[c(1,2)]), x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
	#	} else if(ct.p && !cp.p) {
	#		if(length(ylab)==1) txt <- "Thrust coefficient [-]" else txt <- ylab[2]
	#		if(is.null(leg.text)) leg.text <- c(pc.names, NA, "Power", expression(c[t]))
	#		mtext(txt, 4, line=mgp[1], las=0, cex=cex.lab, col=col.lab)
	#		if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=c(col, NA, rep("darkgray",2)), lty=c(rep(lty[1], length(pc)), NA, lty[c(1,2)]), lwd=c(rep(lwd[1], length(pc)), NA, lwd[c(1,2)]), x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
	#	} else {
	#		if(is.null(leg.text)) leg.text <- c(pc.names)
	#		if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=col, lty=rep(lty[1], length(pc)), lwd=rep(lwd[1], length(pc)), x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
	#	}
	#}
}
