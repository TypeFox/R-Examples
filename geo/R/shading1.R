#' Shading of geoplots?
#' 
#' Shading of geoplots?
#' 
#' 
#' @param cont Contours?
#' @param digits Number of digits?
#' @param colors Colors?
#' @param xlim,ylim Limits?
#' @param fill Fill?
#' @param angle Angle?
#' @param rotate Rotate?
#' @param cex Character expansion?
#' @param rat Ratio?
#' @param minsym Minimum symbol on label?
#' @return No value, shades current geoplot (label?) in some way?
#' @note Needs elaboration.
#' @seealso Called by \code{\link{colsymbol}} and \code{\link{reitaplott}}.
#' @keywords aplot
#' @export shading1
shading1 <-
function(cont, digits, colors, xlim = c(0, 1), ylim = c(0, 1), fill = F, angle,
	rotate, cex, rat, minsym = "<")
{
	xlim <- sort(xlim)
	ylim <- sort(ylim)
	if(cex != 0)
		par(cex = cex)
	ncont <- length(cont)
	if(fill)
		lbox <- max(ncont + 1, 20)
	else lbox <- ncont + 1
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
	xmat[1,  ] <- 0.75
	xmat[2,  ] <- 0.97
	xmat[3,  ] <- 0.97
	xmat[4,  ] <- 0.75
	xmat[5,  ] <- NA
	#       put  text in figure
	par(adj = 0)
	cont <- round(cont, digits = digits)
	textx <- c(1:(length(cont) - 1))
	textx1 <- textx
	textx <- as.character(round(cont[1:(length(cont) - 1)], digits = digits
		))
	textx1 <- as.character(round(cont[2:length(cont)], digits = digits))
	textx <- paste(textx, "-", textx1)
	minsym <- paste(minsym, " ", sep = "")
	textx <- c(paste(minsym, as.character(round(cont[1], digits = digits))),
		textx)
	textx[ncont + 1] <- paste("> ", as.character(round(cont[ncont], digits
		 = digits)))
	boxx <- c(matrix(0.1, 1, length(boxy)))
	boxx <- xlim[1] + (xlim[2] - xlim[1]) * boxx
	boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
	ll <- (ylim[2] - ylim[1]) * 0.05
	if(fill)
		text(boxx, boxy + ll/2, textx)
	else text(boxx, boxy + ll, textx)
	xmat <- xlim[1] + (xlim[2] - xlim[1]) * xmat
	ymat <- ylim[1] + (ylim[2] - ylim[1]) * ymat
	for(i in 1:length(colors)) {
		polygon(xmat[1:4, i], ymat[1:4, i], border = T, density = 
			colors[i], angle = angle)
		angle <- angle + rotate
	}
}

