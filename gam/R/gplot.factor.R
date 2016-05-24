"gplot.factor" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, xlim = NULL, ylim = NULL, fit = TRUE, ...)
{
	if(length(x) != length(y))
		stop("x and y do not have the same length; possibly a consequence of an na.action"
			)
	nn <- as.numeric(table(x))
	codex <- as.numeric(x)
	ucodex <- seq(nn)[nn > 0]
	o <- match(ucodex, codex, 0)
	uy <- as.numeric(y[o])
	ylim <- range(ylim, uy)
	xlim <- range(c(0, sum(nn), xlim))
	rightx <- cumsum(nn)
	leftx <- c(0, rightx[ - length(nn)])
	ux <- ((leftx + rightx)/2)
	delta <- (rightx - leftx)/8
	jx <- runif(length(codex), (ux - delta)[codex], (ux + delta)[codex])
	nnajx <- jx[!is.na(jx)]
	if(rugplot)
		xlim <- range(c(xlim, nnajx))
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
	Levels <- levels(x)
	if(!all(nn>0)) {
		keep <- nn > 0
		ux <- ux[keep]
		delta <- delta[keep]
		leftx <- leftx[keep]
		rightx <- rightx[keep]
		Levels <- Levels[keep]
	}
	plot(ux, uy, ylim = ylim, xlim = xlim, xlab = "", type = "n", ylab = 
		ylab, xaxt = "n", ...)
	mtext(xlab, 1, 2)
	axis(side = 3, at = ux - delta, labels = Levels, srt = 45, tick = FALSE,
		adj = 0)
	if(fit)
		segments(leftx + delta, uy, rightx - delta, uy)
	if(!is.null(residuals))
		points(jx, residuals)
	if(rugplot)
		rug(nnajx)
	if(se) {
		segments(ux + delta, se.upper, ux - delta, se.upper)
		segments(ux + delta, se.lower, ux - delta, se.lower)
		segments(ux, se.lower, ux, se.upper, lty = 2)
	}
	invisible(diff(ylim))
}
