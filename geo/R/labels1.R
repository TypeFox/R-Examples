#' Label plots
#' 
#' Label plots with categories.
#' 
#' 
#' @param cont Contour ?
#' @param digits Number of digits to use in labels
#' @param colors Colors ?
#' @param xlim,ylim Limits ?
#' @param fill Fill with colors?
#' @param minsym Minimum symbol (for lowest category?)?
#' @param label.resolution Label resolution ?
#' @param labtxt Label text ?
#' @param first.color.trans Should first color be transparent? Default TRUE
#' @param mai Margins in inches?
#' @param leftrat Left ratio (giving space for labels??)??
#' @return No value, labels added to current plot.
#' @note Needs elaboration, merge documentation with \code{labels2}, and
#' possibly others?
#' @seealso Called by \code{\link{colsymbol}}, \code{\link{geocontour.fill}}
#' and \code{\link{reitaplott}}.
#' @keywords aplot
#' @export labels1
labels1 <-
function(cont, digits, colors, xlim = c(0, 1), ylim = c(0, 1), fill = F, minsym
	 = "<", label.resolution = 0, labtxt = NULL, first.color.trans = T,
	mai = c(0, 1, 0, 1), leftrat = 0.1)
{
	xlim <- sort(xlim)
	ylim <- sort(ylim)
	dx <- (xlim[2] - xlim[1])
	dy <- (ylim[2] - ylim[1])
	xlim[2] <- xlim[1] + mai[2] * dx
	xlim[1] <- xlim[1] + mai[1] * dx
	ylim[2] <- ylim[1] + mai[4] * dy
	ylim[1] <- ylim[1] + mai[3] * dy
	ncont <- length(cont)
	if(label.resolution == "none")
		lbox <- ncont
	else lbox <- ncont + 1
	if(fill)
		lbox <- max(lbox, 20)
	boxy <- c(1:lbox)
	boxy <-  - boxy/lbox + 1
	boxy1 <- boxy + 1/(1.2 * lbox)
	if(fill) {
		boxy <- boxy[1:(ncont + 1)]
		boxy1 <- boxy1[1:(ncont + 1)]
	}
	ymat <- matrix(0, 5, length(boxy))
	ymat[1,  ] <- boxy
	ymat[2,  ] <- boxy
	ymat[3,  ] <- boxy1
	ymat[4,  ] <- boxy1
	ymat[5,  ] <- NA
	xmat <- matrix(0, 5, length(boxy))
	xmat[1,  ] <- 0.7
	xmat[2,  ] <- 0.95
	xmat[3,  ] <- 0.95
	xmat[4,  ] <- 0.7
	xmat[5,  ] <- NA
	#	put  text in figure
	par(adj = 0)
	cont <- round(cont, digits = digits)
	if(!(label.resolution == "none")) {
		textx <- c(1:(length(cont) - 1))
		textx1 <- textx
		textx <- format(round(cont[1:(length(cont) - 1)] + 
			label.resolution, digits = digits))
		textx1 <- format(round(cont[2:length(cont)], digits = digits))
		textx <- paste(textx, "-", textx1)
		tmp1 <- paste(minsym, format(round(cont[1], digits = digits)))
		tmp2 <- paste(">", format(round(cont[ncont], digits = digits)))
		textx <- c(tmp1, textx, tmp2)
	}
	else {
		print(cont)
		textx <- c(1:length(cont))
		testx <- format(round(cont), digits = digits)
	}
	print(1)
	boxx <- c(matrix(leftrat, 1, length(boxy)))
	boxx <- xlim[1] + abs((xlim[2] - xlim[1])) * boxx
	boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
	ll <- (ylim[2] - ylim[1]) * 0.05
	if(!is.null(labtxt))
		textx <- labtxt
	# put the labels. 
	if(fill) text(boxx, boxy + ll/2, textx) else text(boxx, boxy + ll,
			textx)
	# put the labels. 
	xmat <- xlim[1] + abs((xlim[2] - xlim[1])) * xmat
	ymat <- ylim[1] + (ylim[2] - ylim[1]) * ymat
	if(label.resolution == "none") {
		colors <- colors[2:length(colors)]
	}
	polygon(xmat, ymat, border = T, col = colors)
	if(colors[1] == 0 || first.color.trans) {
		xmat <- c(xmat[1:4], xmat[1])
		# if white color.  
		ymat <- c(ymat[1:4], ymat[1])
		lines(xmat, ymat)
	}
}

