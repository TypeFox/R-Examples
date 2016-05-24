plot.uncertainty <-
function(x, type=c("prob", "uncert"), p.values=c(50, 75, 90), ...) {
###	uncertainty plots
	
	if(missing(x)) stop("uncertainty object 'uncertainty' is mandatory")
	if(missing(type)) type <- "prob"
	type <- match.arg(type)
	
	if(!is.numeric(p.values)) stop("'p.values' must be numeric")
	for(i in 1:length(p.values)) if(p.values[i]%%1!=0 || p.values[i]<1 || p.values[i]>=100) stop("Only positive 'p.values' between 1 and 100 allowed")
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(type=="prob") { # probability of exceedance
		if(any(names(plot.param)=="col")) {
			if(length(plot.param$col)==1) col <- rep(plot.param$col, 4)
			else col <- plot.param$col
		} else {
			if(length(p.values)==1) col <- c("#084081", "#FF0000")
			else if(length(p.values)==2) col <- c("#084081", "#FB6A4A", "#A50F15")
			else if(length(p.values)>2 && length(p.values)<10) {
				if(requireNamespace("RColorBrewer", quietly=TRUE)) col <- c("#084081", RColorBrewer::brewer.pal(length(p.values), "Reds"))
				else col <- c("#084081", rainbow(length(p.values)))
			} else col <- c("#084081", rainbow(length(p.values)))
		}
		if(any(names(plot.param)=="bty")) bty <- plot.param$bty
		else bty <- "o"
		if(any(names(plot.param)=="col.box")) col.box <- plot.param$col.box
		else col.box <- "black"
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
		if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
		else cex.lab <- cex
		if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
		else cex.axis <- cex
		if(any(names(plot.param)=="cex.leg")) cex.leg <- plot.param$cex.leg
		else cex.leg <- cex-0.2
		if(any(names(plot.param)=="lty")) {
			if(length(plot.param$lty)==1) lty <- rep(plot.param$lty, length(p.values)+1)
			else lty <- plot.param$lty
		} else lty <- c(1, rep(2, length(p.values)))
		if(any(names(plot.param)=="lwd")) {
			if(length(plot.param$lwd)==1) lwd <- rep(plot.param$lwd, length(p.values)+1)
			else lwd <- plot.param$lwd
		} else lwd <- c(1.5,rep(1,length(p.values)))
		if(any(names(plot.param)=="x.intersp")) x.intersp <- plot.param$x.intersp
		else x.intersp <- 0.4
		if(any(names(plot.param)=="y.intersp")) y.intersp <- plot.param$y.intersp
		else y.intersp <- 0.8
		if(any(names(plot.param)=="bty.leg")) bty.leg <- plot.param$bty.leg
		else bty.leg <- "n"
		if(any(names(plot.param)=="pos.leg")) pos.leg <- plot.param$pos.leg
		else pos.leg <- "topright"
		if(any(names(plot.param)=="mar")) mar <- plot.param$mar
		else mar <- c(4.5,4.5,1,1)
		if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
		else mgp <- c(2.7,0.7,0)
		if(any(names(plot.param)=="las")) las <- plot.param$las
		else las <- 1
		if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
		else xlab <- "Probability [%]"
		if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
		else ylab <- paste0("Annual energy production [", attr(x$prob.exceedance$aep, "unit"), "]")
		if(any(names(plot.param)=="ylim")) ylim <- plot.param$ylim
		else ylim <- NULL
		if(any(names(plot.param)=="xlim")) xlim <- plot.param$xlim
		else xlim <- NULL
		if(any(names(plot.param)=="legend")) legend <- plot.param$legend
		else legend <- TRUE
		if(any(names(plot.param)=="leg.text")) leg.text <- plot.param$leg.text
		else leg.text <- c("Probability curve", paste("P", p.values, sep=""))
	} else { # uncertainties of methods
		if(any(names(plot.param)=="border")) {
			border <- plot.param$border
			if(length(plot.param$border)==1) border <- rep(plot.param$border, length(x$uncertainty.meth$uncertainty))
			else if(length(plot.param$border)==2) border <- c(rep(plot.param$border[1], length(x$uncertainty.meth$uncertainty)-1), plot.param$border[2])
			else if(length(plot.param$border)==length(x$uncertainty.meth$uncertainty)) border <- plot.param$border
			else stop("Wrong length of border colours")
		} else border <- NA
		if(any(names(plot.param)=="col")) {
			if(length(plot.param$col)==1) col <- rep(plot.param$col, length(x$uncertainty.meth$uncertainty))
			else if(length(plot.param$col)==2) col <- c(rep(plot.param$col[1], length(x$uncertainty.meth$uncertainty)-1), plot.param$col[2])
			else if(length(plot.param$col)==length(x$uncertainty.meth$uncertainty)) col <- plot.param$col
			else stop("Wrong length of colours")
		} else col <- c(rep("#CB181D", length(x$uncertainty.meth$uncertainty)-1), "#940E13")
		if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
		else col.axis <- "black"
		if(any(names(plot.param)=="col.text")) col.text <- plot.param$col.text
		else col.text <- "white"
		if(any(names(plot.param)=="cex")) cex <- plot.param$cex
		else cex <- 1
		if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
		else cex.axis <- cex
		if(any(names(plot.param)=="cex.text")) cex.text <- plot.param$cex.text
		else cex.text <- cex
		if(any(names(plot.param)=="mar")) mar <- plot.param$mar
		else {
			un <- row.names(x$uncertainty.meth)
			l <- nchar(un[1])
			for(n in 1:length(un)) if(nchar(un[n])>l) l <- nchar(un[n])
			mar <- c(1,ceiling(l/2)+(2-round(l/10,1)),1,1)
		}
		if(any(names(plot.param)=="mgp")) mgp <- plot.param$mgp
		else mgp <- c(3,1,0)
		if(any(names(plot.param)=="space")) space <- plot.param$space
		else space <- 0.3
	}
	
	# plot
	if(type=="prob") { # probability of exceedance
		p50 <- attr(x$prob.exceedance$aep, "P50")
		tuc <- tail(x$uncertainty.meth$uncertainty, 1)/100
		p <- qnorm((1-p.values/100), p50, tuc*p50)
		par(mar=mar, mgp=mgp, las=las)
		x <- NULL # just to satisfy R CMD check
		curve(qnorm(rev(x), mean=p50, sd=tuc*p50), xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, col=col[1], lty=lty[1], lwd=lwd[1], cex=cex, cex.lab=cex.lab, xlim=xlim, ylim=ylim, col.axis=col.axis, col.lab=col.lab, bty="n")
		abline(h=p, col=col[2:(length(p.values)+1)], lty=lty[2:(length(p.values)+1)], lwd=lwd[2:(length(p.values)+1)])
		box(bty=bty, col=col.box)
		axis(1, at=seq(0, 1, by=0.2), labels=seq(0, 100, by=20), col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		axis(2, col=col.ticks, col.axis=col.axis, cex.axis=cex.axis)
		if(legend) legend(pos.leg, legend=leg.text, bty=bty.leg, col=col, lty=lty, lwd=lwd, x.intersp=x.intersp, y.intersp=y.intersp, cex=cex.leg, text.col=col.leg)
	} else { # uncertainties of methods
		dat <- rev(x$uncertainty.meth$uncertainty)
		nam <- rev(row.names(x$uncertainty.meth))
		
		par(mar=mar, mgp=mgp, las=1)
		barplot(dat, horiz=TRUE, xaxt="n", yaxt="n", col=rev(col), border=rev(border), space=space)
		bxp <- barplot(dat, horiz=TRUE, space=space, plot=FALSE)
		at <- apply(bxp, 1, mean)
		mtext(nam, side=2, line=mgp[2]-0.5, at=at, cex=cex.axis, col=col.axis)
		text(dat/2, at, paste(dat, "%", sep=""), col=col.text, cex=cex.text)
	}
}
