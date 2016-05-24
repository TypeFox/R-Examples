"gplot.numeric" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, xlim = NULL, ylim = NULL, fit = TRUE, ...)
{
	if(length(x) != length(y))
		stop("x and y do not have the same length; possibly a consequence of an na.action"
			)
### Here we check if its a simple linear term; if so, we include a point at the mean of x
        if(se &&  !is.null(se.y) && ylab==paste("partial for",xlab)){
          x=c(x,mean(x))
          y=c(y,0)
          se.y=c(se.y,0)
                  }
	ux <- unique(sort(x))
	o <- match(ux, x)
	uy <- y[o]
	xlim <- range(xlim, ux)
	ylim <- range(ylim, uy)
	if(rugplot) {
		jx <- jitter(x[!is.na(x)])
		xlim <- range(c(xlim, jx))
	}
	if(se && !is.null(se.y)) {
		se.upper <- uy + 2 * se.y[o]
		se.lower <- uy - 2 * se.y[o]
		ylim <- range(c(ylim, se.upper, se.lower))
	}
	if(!is.null(residuals)) {
		if(length(residuals) == length(y)) {
			residuals <- y + residuals
			ylim <- range(c(ylim, residuals))
		}
		else {
			residuals <- NULL
			warning(paste("Residuals do not match x in \"", ylab,
				"\" preplot object", sep = ""))
		}
	}
	ylim <- ylim.scale(ylim, scale)
	if(!is.null(residuals)) {
		plot(x, residuals, xlim = xlim, ylim = ylim, xlab = xlab, ylab
			 = ylab, ...)
		if(fit)
			lines(ux, uy)
	}
	else {
		if(fit)
			plot(ux, uy, type = "l", xlim = xlim, ylim = ylim,
				xlab = xlab, ylab = ylab, ...)
	}
	if(rugplot)
		rug(jx)
	if(se) {
		lines(ux, se.upper, lty = 3)
		lines(ux, se.lower, lty = 3)
	}
	invisible(diff(ylim))
}
