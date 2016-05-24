#' Fill an area.
#' 
#' The program fills an area. The program is similar to the polygon function
#' except the data is in lat, lon and the transform of the data specified in
#' geoplot is used. Also there are some additional parameters. The graph has to
#' be initialized by geoplot.
#' 
#' 
#' @param lat,lon Latitude and longitude of data ( or x and y coordinates),
#' negative for southern latitudes and western longitudes. May be supplied as
#' two vectors or as a list lat (or x) including vectors lat\$lat and lat\$lon
#' (x\$x and x\$y if projection = none).
#' @param col Color number used.  Default value is 0 (often white).
#' @param border If TRUE borders around the polygon are drawn. Default value is
#' FALSE.
#' @param exterior If TRUE everything that is outside the polygon is painted ,
#' else everything inside. Default value is FALSE. If exterior = TRUE axes and
#' grid often need refreshing by calling geoplot again with new = TRUE.
#' @param nx See geolines for further details.
#' @param outside If TRUE what is outside of the polygon is colored else what
#' is inside. Default value is TRUE.
#' @param plot if TRUE the polygon is plotted. Default is TRUE.
#' @param save if TRUE the points plotted are returned. Default is FALSE.
#' @param rat the ratio of the plot to what is plotted. Default is 0.005
#' meaning that the plot is 0.5\% bigger than what is plotted.
#' @param density see polygon.
#' @param Projection the projection to be used. Default is the one defined by
#' current plot.
#' @param angle see polygon.
#' @param allowed.size printers are limited to printing polygons of certain
#' size, when the polygons size is actually too big for your printer you can
#' enter the tedious task of splitting up the polygon. Default is 4000.
#' @param option Some option. Default 1.
#' @return The points plotted are returned if save = TRUE.
#' @seealso \code{\link{polygon}}, \code{\link{geoplot}},
#' \code{\link{geolines}}, \code{\link{geopoints}}, \code{\link{geotext}},
#' \code{\link{geosymbols}}, \code{\link{geocontour.fill}},
#' \code{\link{geogrid}}, \code{\link{geocontour}}.
#' @examples
#' 
#' \dontrun{     geopolygon(island)              # Paint iceland with 
#'                                      # color #0 (often white).
#' 
#'      geopolygon(island, col = 0, exterior = TRUE)
#' 
#'      geopolygon(geolocator(), col = 1)  # Paints a region defined 
#'                                      # by pointing on map black.
#' 
#'      # Of the maps available island (iceland) is about the only that
#'      # is correctly defined as closed polygon so it is the only one that 
#'      # can be painted by geopolygon.
#' 
#'      geoplot(grid = FALSE, type = "n")
#'      # Star by setting up the plot.
#'      geopolygon(gbdypif.500, col = 4, exterior = FALSE, r = 0)
#'      # Use geopolygon to draw the 500 m area. 
#'      geopolygon(gbdypif.100, col = 155, exterior = FALSE, r = 0)
#'      # Draw 100 m are over the 500 m. 
#'      geolines(eyjar, col = 115)
#'      # Add islands around Iceland.
#'      gbplot(c(100, 500), depthlab = TRUE)
#'      # Draw the depth lines, labels on lines.
#'      geopolygon(island, col = 115, outside = TRUE, r = 0)
#'      # Draw Iceland over.
#'      geoplot(grid = FALSE, new = TRUE)
#'      # Draw lines around Iceland, could also use geolines.
#' }
#' @export geopolygon
geopolygon <-
function(lat, lon = NULL, col = "white", border = FALSE, exterior = FALSE, nx = 1,
	outside = FALSE, plot = TRUE, save = FALSE, rat = 0.005, density = -1, Projection
	 = NULL, angle = 45, allowed.size = 80000, option = 1)
{
	geopar <- getOption("geopar")
	if(is.null(Projection))
		Projection <- geopar$projection
	# 	for structures too large for hardware
	index <- lat$index
	RANGE <- lat$range
	LENGTH <- lat$length
	if(exterior)
		in.or.out <- 1
	else in.or.out <- 0
	err <- FALSE
	if(is.null(lon)) {
		if(Projection == "none") {
			lon <- lat$y
			lat <- lat$x
		}
		else {
			lon <- lat$lon
			lat <- lat$lat
		}
	}
	if(Projection != "none") {
		# degrees and minutes
		if(mean(lat, na.rm = TRUE) > 1000) {
			lat <- geoconvert(lat)
			lon <-  - geoconvert(lon)
		}
	}
	if(length(lat) == 2) {
		lat <- c(lat[1], lat[1], lat[2], lat[2], lat[1])
		lon <- c(lon[1], lon[2], lon[2], lon[1], lon[1])
	}
	if(nx > 1) {
		# fill in with points for lambert. 
		x <- fill.points(lat, lon, nx, option = 2)
		lat <- x$x
		lon <- x$y
	}
	oldpar <- selectedpar()
	par(geopar$gpar)
	if(outside)
	par(xpd = TRUE)
	else par(xpd = FALSE)
	on.exit(par(oldpar))
	gx <- geopar$limx
	rx <- gx[2] - gx[1]
	gy <- geopar$limy
	ry <- gy[2] - gy[1]
	gx[1] <- gx[1] + rat * rx
	gx[2] <- gx[2] - rat * ry
	gy[1] <- gy[1] + rat * ry
	gy[2] <- gy[2] - rat * ry
	brd <- data.frame(x = c(gx[1], gx[2], gx[2], gx[1], gx[1]), y = c(
		gy[1], gy[1], gy[2], gy[2], gy[1]))
	brd1 <- invProj(brd)
	brd1 <- data.frame(lat = brd1$lat, lon = brd1$lon)
	if(!is.null(index)) {
		limits <- invProj(geopar$limx, geopar$limy)
		for(i in 1:length(index)) {
			xx <- Proj(lat[index[[i]]], lon[index[[i]]])
			if(!outside)
				xx <- cut_multipoly(xx, brd, in.or.out)
			if(length(xx$x) > 0)
				polygon(xx$x, xx$y, col = col, border = border,
					density = density, angle = angle)
		}
		return(invisible())
	}
	else {
		if(exterior) {
			i1 <- geoinside(brd1, data.frame(lat = lat, lon = lon),
				option = 0)
			i <- 1:4
			i <- i[is.na(match(i, i1))]
			if(length(i) == 0)
				return(invisible())
			else i1 <- i[1]
			i <- geoinside(data.frame(lat = lat, lon = lon), brd1,
				na.rm = TRUE, robust = FALSE, option = 0)
			if(length(i) == length(lat) || option != 1) {
				lat <- lat[!is.na(lat)]
				lon <- lon[!is.na(lon)]
				dist <- (lat - brd1$lat[i1])^2 + (lon - brd1$
					lon[i1])^2 * cos((mean(lat) * pi)/
					180)^2
				o <- order(dist)
				lat <- c(lat[c(o[1]:length(lat), 1:o[1])],
					brd1$lat[c(i1:4, 1:i1)], lat[o[1]])
				lon <- c(lon[c(o[1]:length(lon), 1:o[1])],
					brd1$lon[c(i1:4, 1:i1)], lon[o[1]])
				xx <- Proj(lat, lon, geopar$scale, geopar$
					b0, geopar$b1, geopar$l1, Projection)
				if(plot) {
					polygon(xx$x, xx$y, col = col, border
						 = border, density = density,
						angle = angle)
					return(invisible())
				}
				else return(invisible(invProj(xx)))
			}
		}
	}
	err <- FALSE
	xx <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		Projection)
	if(!outside)
		xx <- cut_multipoly(xx, brd, in.or.out)
	if(length(xx$x) > allowed.size && plot) {
		ind <- seq(along = xx$x)
		ind <- ind[is.na(xx$x)]
		if(length(ind) == 0)
			err <- TRUE
		else {
			ind <- c(1, ind, length(xx$x))
			if(max(diff(ind)) > allowed.size)
				err <- TRUE
			else err <- FALSE
		}
	}
	if(plot) {
		if(!err)
			polygon(xx$x, xx$y, col = col, border = border, density
				 = density, angle = angle)
		else print("too large polygon, change parameter allowed.zize")
	}
	if(save) {
		xx <- invProj(xx$x, xx$y, geopar$scale, geopar$b0, geopar$
			b1, geopar$l1, Projection)
		return(list(lat = xx$lat, lon = xx$lon))
	}
	return(invisible())
}

