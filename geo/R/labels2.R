#' Label plots
#' 
#' Label plots with categories.
#' 
#' 
#' @param cont Contour ?
#' @param digits Number of digits to use in labels
#' @param colors Colors ?
#' @param xlim,ylim Limits ?
#' @param nx ??. Default 4
#' @param fill Fill ?
#' @return No value, labels added to current plot.
#' @note Needs elaboration, merge documentation with \code{labels1}, and
#' possibly others?
#' @seealso alled by \code{\link{colsymbol}} and \code{\link{geocontour.fill}}.
#' @keywords aplot
#' @export labels2
labels2 <-
function(cont, digits, colors, xlim = c(0, 1), ylim = c(0, 1), nx = 4, fill = F
	)
{
	xlim <- sort(xlim)
	ylim <- sort(ylim)
	ncont <- length(cont)
	lbox <- ncont + 1
	if(fill)
		lbox <- max(lbox, 20)
	boxy <- c(1:lbox)
	boxy <-  - boxy/(lbox + 2) + 1
	dy <- 1/lbox
	boxy <- boxy - dy/2
	boxy1 <- boxy + 1/lbox
	ymat <- matrix(0, 5, length(boxy))
	ymat[1,  ] <- boxy
	ymat[2,  ] <- boxy
	ymat[3,  ] <- boxy1
	ymat[4,  ] <- boxy1
	ymat[5,  ] <- NA
	xmat <- matrix(0, 5, length(boxy))
	xmat[1,  ] <- 0.6
	xmat[2,  ] <- 0.9
	xmat[3,  ] <- 0.9
	xmat[4,  ] <- 0.6
	xmat[5,  ] <- NA
	#	put  text in figure
	ind <- c(1, c(1:floor((length(cont))/nx)) * nx)
	if(ind[length(ind)] == (length(cont)))
		ind <- c(ind, (length(cont)))
	par(adj = 0)
	cont <- round(cont, digits = digits)
	textx <- format(round(cont[ind], digits = digits))
	boxx <- c(matrix(0.1, 1, length(boxy)))
	boxx <- xlim[1] + (xlim[2] - xlim[1]) * boxx
	boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
	text(boxx[ind], boxy[ind], textx)
	# put the lables.  
	xmat <- xlim[1] + abs((xlim[2] - xlim[1])) * xmat
	ymat <- ylim[1] + (ylim[2] - ylim[1]) * ymat
	polygon(xmat, ymat, border = F, col = colors)
	if(colors[1] == 0) {
		xmat <- c(xmat[1:4], xmat[1])
		# if white color.  
		ymat <- c(ymat[1:4], ymat[1])
		lines(xmat, ymat)
	}
}

