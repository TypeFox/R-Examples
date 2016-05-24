`drawkpc` <-
function(model, err = TRUE, pointcol = "blue", rdcol = "red", noisecol = "black", errcol = "brown", ...)
{
	xlab <- "kernel pca coefficients"
	ylab <- ""
	x <- 1:length(model$kpc)
	y <- abs(model$kpc)
	legend <- "relevant dimension"
	col <- rdcol

	plot(x = x, y = y, xlab = xlab, ylab = ylab, col = pointcol, ...)
	
	# draw noise, if available
	if(!is.null(model$noise))
	{
		segments(x0 = 0, y0 = model$noise, x1 = max(x), y1 = model$noise, col = noisecol, ...)
		legend <- c(legend, "noise")
		col <- c(col, noisecol)
	}
	
	# draw loo-cv-error/loglik?
	if(err)
	{
		errheight <- max(model$err) - min(model$err)
		kpcheight <- max(y) - min(y)
		scalederr <- 0.8*kpcheight*(model$err/errheight)
		align <- max(y) - max(scalederr)
		lines(x = 1:length(scalederr), y = scalederr + align, col = errcol)
		if(model$tcm)
		{
			legend <- c(legend, "negative loglik (scaled)")
		}
		else
		{
			legend <- c(legend, "loo-cv-error (scaled)")
		}
		col <- c(col, errcol)
	}
	
	# draw relevant dimension
	segments(x0 = model$rd, y0 = 0, x1 = model$rd, y1 = max(y), col = rdcol, ...)
	
	# legend
	legend(x = length(model$kpc)/2, y = max(abs(model$kpc)), legend = legend, col = col, lwd = 2, bg = "white")
}

