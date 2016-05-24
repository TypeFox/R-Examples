#' Labels line ?
#' 
#' Labels line ?.
#' 
#' 
#' @param cont Contours?
#' @param digits Number of digits
#' @param colors Colors
#' @param lty Line types
#' @param lwd Line widths
#' @param xlim,ylim Limit
#' @param linew Linewidth for some reason not lwd?
#' @return No value, labels added to current plot.
#' @note Needs elaboration, possibly merge documentation with others label
#' functions?
#' @seealso Called by \code{\link{geocontour}}.
#' @keywords aplot
#' @export labels_line
labels_line <-
function(cont, digits, colors, lty, lwd, xlim = c(0, 1), ylim = c(0, 1), linew = F)
{
	xlim <- sort(xlim)
	ylim <- sort(ylim)
	ncont <- length(cont)
	if(length(lty) == ncont)
		linetypes <- T
	else linetypes <- F
	lbox <- ncont
	boxy <- c(1:lbox)
	boxy <-  - boxy/(lbox + 1) + 1
	boxy1 <- boxy + 1/(1.2 * lbox)
	ymat <- matrix(0, 2, length(boxy))
	ymat[1,  ] <- boxy
	ymat[2,  ] <- boxy
	xmat <- matrix(0, 2, length(boxy))
	xmat[1,  ] <- 0.7
	xmat[2,  ] <- 0.95
	#	put  text in figure
	par(adj = 0)
	cont <- round(cont, digits = digits)
	textx <- format(cont)
	boxx <- c(matrix(0.1, 1, length(boxy)))
	boxx <- xlim[1] + abs((xlim[2] - xlim[1])) * boxx
	boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
	ll <- (ylim[2] - ylim[1]) * 0.04
	text(boxx, boxy + ll, textx, col = 1)
	# put the lables.  
	xmat <- xlim[1] + abs((xlim[2] - xlim[1])) * xmat
	ymat <- ylim[1] + (ylim[2] - ylim[1]) * ymat
	for(i in 1:ncont) {
		if(linew)
			par(lwd = lwd[i])
		if(linetypes)
			par(lty = lty[i])
		lines(xmat[, i], ymat[, i] + ll, col = colors[i])
	}
}

