#' Performs the inverse Mercator or Lambert projection of data.
#' 
#' Accepts data in spherical coordinates with center lattitude 0 and longitude
#' 0 and performs an inverse Mercator, Lambert or no tranformation.
#' 
#' 
#' @param x,y The input data to be inversely projected, may be given as two
#' vectors or as list attributes (\$x and \$y).
#' @param scale The scale of the input date (m, km or miles), default is the
#' scale defined in geopar (the scale defined when the plot is initialized).
#' @param b0 if projection = Mercator b0 is the center of the Mercator
#' projection. If projection = Lambert b0, b1 are the latitudes defining the
#' Lambert projection. Default are the b0 and b1 defined in geopar.
#' @param b1 Second defining latitude for Lambert projection.
#' @param l1 The longitude defining the Lambert projection, default is the l1
#' defined in geopar.
#' @param projection The projection to be inversed, legal projections are
#' "mercator", "Lambert" and "none". Default is the projection defined in
#' geopar.
#' @return The function returns a list containing if projection = "none" \$x
#' and \$y, if projection is mercator or Lambert it includes the projection
#' (\$projection), the scale (\$scale), \$lat and \$lon and \$x and \$y.
#' @seealso \code{\link{invProj}}, \code{\link{geopar}}, \code{\link{geoplot}}.
#' @export invProj
invProj <-
function(x, y = NULL, scale = getOption("geopar")$scale,
         b0 = getOption("geopar")$b0, b1 = getOption("geopar")$b1,
         l1 = getOption("geopar")$l1,
         projection = getOption("geopar")$projection)
{
	if(is.null(y)) {
		y <- x$y
		x <- x$x
	}
	if(projection == "Lambert") {
		x <- invlambert(x, y, b0, l1, b1, scale, old = T)
	}
	else if(projection == "Mercator") {
		x <- invmerc(x, y, scale, b0)
	}
	else if(projection == "none") {
		x <- list(x = x, y = y)
	}
}

