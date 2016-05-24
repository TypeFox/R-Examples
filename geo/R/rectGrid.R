#' Produce a grid of rectangles on a plot, filled with colors if desired.
#' 
#' 
#' Functions that plot a line grid or one filled with colors for rectangle
#' codes in systems of various resolutions with call to \code{\link{geolines}}
#' or \code{\link{geopolygon}} and perimeter utility functions (see
#' \code{\link{rectPeri}}).
#' 
#' The functions simply outline or fill the rectangles they are given.  whereas
#' \code{\link{reitaplott}} and \code{\link{geoSR}} assume levels and are more
#' hi-level.
#' 
#' @name rGrid
#' @aliases rgrid srgrid mrgrid drgrid
#' @param r,sr,dr,mr Codes of rectangle to be outlined or filled. In different
#' resolutions, \code{r, sr} with their own system (see
#' \code{\link{deg2rect}}), \code{mr} dimensions based on minutes, \code{dr} on
#' degrees.
#' @param dlat,dlon Dimensions of latitude and longitude given in minutes and
#' degrees for \code{mrgrid} and \code{drgrid}, respectively
#' @param fill Logical, whether or not to fill the plotted rectangles.
#' @param \dots other arguments to \code{\link{geopolygon}} or
#' \code{\link{geolines}} as appropriate, notably \code{col}.
#' @return No values returned, used for side-effects.
#' @note Functions \code{\link{reitaplott}} and \code{\link{geoSR}} have more
#' in-built functionality to deal with level-plots of rectangles.
#' @author STJ
#' @seealso \code{\link{deg2rect}}, \code{\link{geolines}},
#' \code{\link{geopolygon}}, \code{\link{rectPeri}}, \code{\link{reitaplott}},
#' \code{\link{geoSR}}.
#' @keywords hplot spatial
#' @examples
#' 
#' 
#' geoplot(grid = FALSE)
#' tmp <- island
#' tmp$sr <- d2sr(island) 
#' srects <- aggregate(. ~ sr, tmp, length)
#' names(srects)[2] <- "count"
#' srects$lev <- cut(srects$count, c(0, 1, 5, 10, 20, 50, 100))
#' mycol <- heat.colors(length(unique(srects$lev)))
#' srgrid(srects$sr, fill = TRUE, col = mycol[srects$lev])
#' geolines(island)
#' 
#'

#' @export rgrid
#' @rdname rectGrid
rgrid <-
function(r, fill = FALSE, ...)
{
	n <- length(r)
	lat <- r2d(r)$lat
	lon <- r2d(r)$lon
	lat <- c(rep(lat, 5), rep(NA, n))
	lon <- c(rep(lon, 5), rep(NA, n))
	lat <- as.vector(matrix(matrix(lat, nrow = 6, byrow = T), ncol = 1))
	lon <- as.vector(matrix(matrix(lon, nrow = 6, byrow = T), ncol = 1))
	lat <- lat - c(-1/4, 1/4, 1/4, -1/4, -1/4, NA)
	lon <- lon - c(-0.5, -0.5, 0.5, 0.5, -0.5, NA)
	lat <- data.frame(lat = lat, lon = lon)
	if(fill)
		geopolygon(lat, ...)
	else geolines(lat, ...)
}

#' @export srgrid
#' @rdname rectGrid
srgrid <-
function(sr, fill = FALSE, ...)
{
	n <- length(sr)
	lat <- sr2d(sr)$lat
	lon <- sr2d(sr)$lon
	lat <- c(rep(lat, 5), rep(NA, n))
	lon <- c(rep(lon, 5), rep(NA, n))
	lat <- as.vector(matrix(matrix(lat, nrow = 6, byrow = T), ncol = 1))
	lon <- as.vector(matrix(matrix(lon, nrow = 6, byrow = T), ncol = 1))
	lat <- lat - c(-1/8, 1/8, 1/8, -1/8, -1/8, NA)
	lon <- lon - c(-0.25, -0.25, 0.25, 0.25, -0.25, NA)
	lat <- data.frame(lat = lat, lon = lon)
	if(fill)
		geopolygon(lat, ...)
	else geolines(lat, ...)
}

#' @export mrgrid
#' @rdname rectGrid
mrgrid <-
function(mr, dlat = 5, dlon = 10, fill = FALSE, ...)
{
  n <- length(mr)
  lat <- mr2d(mr, dlat = dlat, dlon = dlon)$lat
  lon <- mr2d(mr, dlat = dlat, dlon = dlon)$lon
  lat <- c(rep(lat, 5), rep(NA, n))
  lon <- c(rep(lon, 5), rep(NA, n))
  lat <- as.vector(matrix(matrix(lat, nrow = 6, byrow = T), ncol = 1))
  lon <- as.vector(matrix(matrix(lon, nrow = 6, byrow = T), ncol = 1))
  lat <- lat + c(dlat/120, dlat/120,  - dlat/120,  - dlat/120, dlat/
  	120, NA)
  lon <- lon + c(dlon/120,  - dlon/120,  - dlon/120, dlon/120, dlon/
  	120, NA)
  lat <- data.frame(lat = lat, lon = lon)
  if(fill)
  	geopolygon(lat, ...)
  else geolines(lat, ...)
}

#' @export drgrid
#' @rdname rectGrid
drgrid <-
function(dr, dlat = 1, dlon = 2, fill = FALSE, ...)
{
  n <- length(dr)
  lat <- dr2d(dr, dlat = dlat, dlon = dlon)$lat
  lon <- dr2d(dr, dlat = dlat, dlon = dlon)$lon
  lat <- c(rep(lat, 5), rep(NA, n))
  lon <- c(rep(lon, 5), rep(NA, n))
  lat <- as.vector(matrix(matrix(lat, nrow = 6, byrow = T), ncol = 1))
  lon <- as.vector(matrix(matrix(lon, nrow = 6, byrow = T), ncol = 1))
  lat <- lat + c(dlat/2, dlat/2,  - dlat/2,  - dlat/2, dlat/
  	2, NA)
  lon <- lon + c(dlon/2,  - dlon/2,  - dlon/2, dlon/2, dlon/
  	2, NA)
  lat <- data.frame(lat = lat, lon = lon)
  if(fill)
  	geopolygon(lat, ...)
  else geolines(lat, ...)
}

