#' Define regions intteractively
#' 
#' Define regions interactively.
#' 
#' 
#' @param nholes The number of holes in the data, number of regions - 1
#' @return Returns a composite list of the following compents: \item{reg}{List
#' with components number of components (?)} \item{x, y}{Coordinates when
#' \code{geopar$projection == "none"}} \item{lat, lon}{Geographical coordinates
#' for other projections (Mercator and Lambert (?)} \item{lxv}{ID of region
#' (?)}
#' @note Needs further elaboration.
#' @seealso Called by \code{\link{geodefine}}, calls \code{\link{invProj}}.
#' @keywords iplot
#' @export giveborder
giveborder <-
function(nholes = 0)
{
	geopar <- getOption("geopar")
	cat(" Give the border of estimating area (not holes) : \n")
	xb <- locator(type = "l")
	lines(c(xb$x[length(xb$x)], xb$x[1]), c(xb$y[length(xb$x)], xb$y[1]))
	# last to fyrst point
	xba <- xb
	# vector to use for plotting
	xba$x <- c(xb$x, xb$x[1])
	#new
	xba$y <- c(xb$y, xb$y[1])
	#new
	#	Read in the data for holes.  Stored in two vectors (added to
	#	the border.  xba has NA's and the first point of each point
	#	twice and is used for plotting while xb is used for the
	#	program "margh" and has each point only once and no NA's.
	lxv <- c(1:(nholes + 2))
	lxv[] <- 0
	if(nholes > 0) {
		lx1 <- 1
		#index
		for(i in 1:nholes) {
			xba$x <- c(xba$x, xb$x[lx1], NA)
			# vector to use in
			xba$y <- c(xba$y, xb$y[lx1], NA)
			# contourplots first point	    
			cat("Hole (or other area) ", i, "\n")
			lx1 <- length(xb$x) + 1
			# preserve last length
			lxv[i + 1] <- length(xb$x)
			xb1 <- locator(type = "l")
			xb$x <- c(xb$x, xb1$x)
			xb$y <- c(xb$y, xb1$y)
			xba$x <- c(xba$x, xb1$x)
			xba$y <- c(xba$y, xb1$y)
			lx2 <- length(xb$x)
			xba$x <- c(xba$x, xb$x[lx1])
			xba$y <- c(xba$y, xb$y[lx1])
			#new
			lines(c(xb$x[lx2], xb$x[lx1]), c(xb$y[lx2], xb$y[lx1]))
		}
	}
	if(nholes == 0)
		lx1 <- 1
	#  	xba$x <-c(xba$x,xb$x[lx1])
	#	xba$y <-c(xba$y,xb$y[lx1])
	reg <- invProj(xba$x, xba$y, geopar$scale, geopar$b0, geopar$b1, geopar$
		l1, geopar$projection)
	lxv[nholes + 2] <- length(xb$x)
	xb <- invProj(xb$x, xb$y, geopar$scale, geopar$b0, geopar$b1, geopar$
		l1, geopar$projection)
	if(geopar$projection == "none") {
		x <- xb$x
		y <- xb$y
		return(list(reg = reg, x = x, y = y, lxv = lxv))
	}
	else {
		lat <- xb$lat
		lon <- xb$lon
		return(list(reg = reg, lat = lat, lon = lon, lxv = lxv))
	}
}

