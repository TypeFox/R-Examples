#' Produce a grid over an area
#' 
#' Produce a grid over an area, possibly interactively.
#' 
#' 
#' @param lat Latitude
#' @param lon Longitude, if not included in \code{lat}
#' @param type Plot type
#' @param pch Plot character
#' @param xlim,ylim Limits for plot
#' @param b0 Base latitude
#' @param r Range expansion
#' @param country Country plotted
#' @param xlab,ylab Labels for x and y axes
#' @param option Method for determining plotted area extent, default "cut" (to
#' the range of the data)
#' @param reg Region to be gridded, can be set interactively
#' @param dx Resolution in each direction (?)
#' @param nx Number of gridpoints in each direction (?)
#' @param grpkt Gridpoints can also be supplied to the function for plotting
#' (??????????)
#' @param scale Projection scale (general \code{geo} default is "km")
#' @param find Should the gridpoints within \code{reg} be determined?
#' @param new Plot control argument ?
#' @param grid Draw grid (which grid?)
#' @param projection Projection to use
#' @param n Number of gridpoints (?)
#' @param b1 Second latitude for the Lambert projection
#' @param nholes number of holes to be sent to \code{geodefine} when setting
#' out the region to be gridded
#' @return List of components: \item{grpt}{Gridpoints} \item{reg}{Region over
#' which the grid was laid} \item{find}{Was \code{find = TRUE}?} \item{xgr}{If
#' \code{find = TRUE} the gridpoints within region \code{reg} is also
#' returned.}
#' @note Needs elaboration, check use of \code{find = TRUE}
#' @seealso Calls a number of functions, i.e. \code{\link{geodefine}},
#' \code{\link{geoplot}}, \code{\link{geopoints}}, \code{\link{gridpoints}},
#' \code{\link{inside}}, \code{\link{selectedpar}}
#' @keywords aplot iplot
#' @export setgrid
setgrid <-
function(lat, lon = 0, type = "p", pch = "*", xlim = c(0, 0), ylim = c(0, 0),
	b0 = 65, r = 1.1, country = geo::island, xlab = "default", ylab = "default",
	option = "cut", reg = 0, dx = c(0, 0), nx = c(0, 0), grpkt = 0, scale
	 = "km", find = F, new = F, grid = T, projection = "Mercator", n = 2500,
	b1 = b0, nholes = 0)
{
	geopar <- getOption("geopar")
	if(length(lon) == 1) {
		if(projection == "none") {
			lon <- lat$y
			lat <- lat$x
		}
		else {
			lon <- lat$lon
			lat <- lat$lat
		}
	}
	geoplot(lat, lon, type = type, pch = pch, xlim = xlim, ylim = ylim,
		b0 = b0, r = r, country = country, xlab = xlab, ylab = ylab,
		option = option, new = new, grid = grid, projection = 
		projection, b1 = b1)
	# Find borders either given or with the locator.  
	oldpar <- selectedpar()
	par(geopar$gpar)
	# set graphical parameters
	on.exit(par(oldpar))
	if(length(reg) == 1) {
		# use the locator.  
		reg <- geodefine(nholes = nholes)
	}
	xgr <- gridpoints(reg, dx, grpkt, nx, n)
	# grid points.  
	grpt <- xgr$xgr
	xgr <- xgr$xgra
	# change names
	geopoints(xgr, pch = ".")
	# 	Find what is inside the borders.  
	if(find) {
		xgr <- inside(xgr, reg = reg)
		geopoints(xgr, pch = "+")
		return(list(xgr = xgr, grpt = grpt, reg = reg, find = find))
	}
	else return(list(grpt = grpt, reg = reg, find = find))
}

