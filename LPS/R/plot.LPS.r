## Linear Predictor Score plot
## Author : Sylvain Mareschal <maressyl@gmail.com>
plot.LPS <- function(
		x,
		y,
		method = c("Wright", "Radmacher", "exact"),
		threshold = 0.9,
		values = FALSE,
		col.classes = c("#FFCC00", "#1144CC"),
		xlim,
		yaxt = "s",
		xlab = "LPS",
		ylab,
		las = 0,
		lwd = 2,
		...
	) {
	# X
	if(missing(xlim)) xlim <- c(min(x$means - 4*x$sds), max(x$means + 4*x$sds))
	xval <- seq(from=xlim[1], to=xlim[2], length.out=1000)
	
	# Y
	method <- match.arg(method)
	if(missing(y)) y <- "density"
	y <- match.arg(y, c("density", "probability"))
	if(missing(ylab)) ylab <- sprintf("%s (%s)", y, method)
	if(y == "density") {
		# Density plot
		if(method == "Wright") {
			# Wright et al (gaussian densities)
			yval1 <- dnorm(xval, mean=x$means[1], sd=x$sds[1])
			yval2 <- dnorm(xval, mean=x$means[2], sd=x$sds[2])
			ylim <- range(c(yval1, yval2), na.rm=TRUE)
			type <- "l"
		} else {
			# Other (discrete values)
			xval <- c(x$scores[[1]], x$scores[[2]])
			yval1 <- c(rank(x$scores[[1]]), rep(NA, length(x$scores[[2]])))
			yval2 <- c(rep(NA, length(x$scores[[1]])), rank(x$scores[[2]]))
			ylim <- range(c(yval1, yval2), na.rm=TRUE)
			type <- "p"
		}
	} else {
		# Probability plot
		p <- predict(x, newdata=xval, type="probability", method=method, plot=FALSE)
		yval1 <- p[,1]
		yval2 <- p[,2]
		ylim <- 0:1
		type <- "l"
	}
	
	# Compute Y ticks on true ylim
	if(yaxt == "s") yat <- pretty(ylim)
	
	# Add Y space to plot values
	if(isTRUE(values)) {
		ypoints <- ylim[1] - (ylim[2] - ylim[1]) * c(0.1, 0.2)
		ylim[1] <- ylim[1] - (ylim[2] - ylim[1]) * 0.25
	}
	
	# Plot (yaxt handled manually)
	plot(x=xval, y=yval1, type=type, col=col.classes[1], xlim=xlim, ylim=ylim, yaxt="n", xlab=xlab, ylab=ylab, las=las, lwd=lwd, ...)
	par(new=TRUE)
	plot(x=xval, y=yval2, type=type, col=col.classes[2], xlim=xlim, ylim=ylim, yaxt="n", xlab="",   ylab="",   las=las, lwd=lwd, ...)
	
	# Radmacher's means
	if(method == "Radmacher") {
		abline(v=x$means[1], col=col.classes[1])
		abline(v=x$means[2], col=col.classes[2])
	}
	
	# Gray zone (score thresholds)
	if(!is.na(threshold) & method != "Radmacher") {
		if(y == "density") p <- predict(x, newdata=xval, type="probability", method=method, plot=FALSE)
		if(p[1,1] < p[nrow(p),1]) {
			grayzone <- c(
				mean(xval[  0:1 + tail(which(p[,1] < threshold), 1) ]),
				mean(xval[ -1:0 + head(which(p[,2] < threshold), 1) ])
			)
		} else {
			grayzone <- c(
				mean(xval[  0:1 + tail(which(p[,2] < threshold), 1) ]),
				mean(xval[ -1:0 + head(which(p[,1] < threshold), 1) ])
			)
		}
		rect(xleft=min(grayzone), xright=max(grayzone), ybottom=-1, ytop=1, col="#00000033", border=NA)
		mtext(text=sprintf("%.3f", grayzone), at=grayzone, side=1)
		mtext(text=sprintf("p < %g%%", threshold*100), at=mean(grayzone), side=3)
	}
	
	# Values
	if(isTRUE(values)) {
		points(x=x$scores[[1]], y=rep(ypoints[1], length(x$scores[[1]])), col=col.classes[1], pch="|")
		points(x=x$scores[[2]], y=rep(ypoints[2], length(x$scores[[2]])), col=col.classes[2], pch="|")
	}
	
	# Y axis
	if(yaxt == "s") axis(side=2, at=yat, labels=yat, las=las)
	
	# Legend
	legend(
		x = "topleft",
		inset = 0.01,
		bg = "#EEEEEE",
		legend = c(x$classes),
		col = col.classes,
		lty = "solid",
		lwd = lwd
	)
}

