#' Performs Mercator or Lambert projection of data.
#' 
#' Performs Mercator, Lambert or no projection of data, as a default it will
#' perform the projection used in current plot.
#' 
#' 
#' @param a,b The input data to be projected, may be given as two vectors or as
#' list attributes, \$lat and \$lon (\$x and \$y if projection = none).
#' @param scale The scale used for the projection, (m, km or miles). Default is
#' the scale defined in geopar (the scale defined when the plot is
#' initialized).
#' @param b0 if projection = mercator b0 is the center of the mercator
#' projection. If projection = "Lambert" b0 and b1 are lattitudes defining the
#' Lambert projection. Default are the b0 and b1 defined in geopar.
#' @param b1 Second defining latitute for Lambert projection
#' @param l1 The longitude defining the Lambert projection, default is the l1
#' defined in geopar.
#' @param projection The projection of the data, legal projections are
#' "Mercator", "Lambert" and "none".
#' @param col.names This has to be set to the default value of c("lon", "lat"),
#' otherwise projection will be set to "none".
#' @return The function returns a list containing if projection = "none" x and
#' y, if projection is "Mercator" or "Lambert" it includes the projection
#' (\$projection), the scale (\$scale), \$lat and \$lon and \$x and \$y (the
#' distance in \$scale from point (0,0) in spherical coordinates.
#' @seealso \code{\link{invProj}}, \code{\link{geopar}}, \code{\link{geoplot}}.
#' @examples
#' 
#'   # For an example of use for this function see i.e. init() where
#'   # it is called:
#' \dontrun{
#'   xgr <- Proj(lat, lon, scale, b0, b1, l1, projection)
#' }
#' 
#' @export Proj
Proj <-
function(a, b = 0, scale = getOption("geopar")$scale,
         b0 = getOption("geopar")$b0, b1 = getOption("geopar")$b1,
         l1 = getOption("geopar")$l1,
         projection = getOption("geopar")$projection,
         col.names = c("lon", "lat"))
{
	if(col.names[1] != "lon" || col.names[2] != "lat")
		projection <- "none"
	if(is.list(a)) {
		if(projection == "none") {
			b <- a$y
			a <- a$x
		}
		else {
			b <- a$lon
			a <- a$lat
		}
	}
	if(projection == "Lambert") {
		x <- lambert(a, b, b0, l1, b1, scale, old = T)
	}
	else if(projection == "Mercator") {
		x <- mercator(a, b, scale, b0)
	}
	else if(projection == "none") {
		x <- list(x = a, y = b)
	}
}

