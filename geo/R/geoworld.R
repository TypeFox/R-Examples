#' Plots rough outline of the world.
#' 
#' Plots roughly the outlines of those countries which have borders
#' intersecting the current plot initialized by geoplot.  The countries can
#' also be filled.
#' 
#' 
#' @param regions Allows plotting only a certain part of the world to speed up
#' the function, such as regions = "iceland". Default is plotting what is
#' insidie the map.
#' @param exact Draws a more exact plot.
#' @param boundary If true country boundaries are drawn.  Default is false.
#' @param fill If true countries are filled with a color.  Default is false.
#' @param color Default is 1 (usually black).
#' @param lwd Linewidth, default is 1.
#' @param lty Linetype, default is 1.
#' @param plot If false a plot is not drawn.  Default is true.
#' @param type If lines or points should be potted, see type for geoplot.
#' @param pch The character to be used for plotting, see par.
#' @param database The database the world plot is to be taken from.  Currently
#' there are two possible databases world and worldHires.  world is in the map
#' library and worlHires in the mapdata library and is much more precies.
#' @param return.data If true those points used to make the in the plot are
#' return, default is FALSE. Can be used to save the coastlines of specified
#' countries.
#' @param outside If TRUE data are plotted outside the borders.  Default is
#' FALSE.
#' @param allowed.size Maximum size of polygon that can be filled.  Default is
#' 80 000 which is rather large polygon but this value is rapidly changing by
#' more powerful hardware.
#' @return default none, see option return data.
#' @section Side Effects: the outlines of all countries who intersect current
#' plot are drawn.
#' @seealso \code{\link{geoplot}}, \code{\link{fill.outside.border}},
#' \code{\link{geopolygon}}, \code{\link{par}}.
#' @examples
#' 
#' \dontrun{     geoplot(xlim = c(0, -53), ylim = c(53, 75))
#'      geoworld()
#' 
#'      # Should plot in all countries who intersect the plot draw
#'      # with geoplot.
#' 
#' # The packages maps and mapdata need to be installed
#' # worldHires is a very detailed database of coastlines from the 
#' # package mapdata.  Could be problematic if used with fill = TRUE)
#' # Allowed.size is the maximum allowed size of polygons.  
#' library(map) # world coastlines and programs 
#' library(mapdata) # more detailed coastlines
#' geoplot(xlim = c(20, 70), ylim = c(15, 34))
#' geoworld(database = "worldHires", fill = TRUE, col = 30, allowed.size = 1e5)
#' 
#' geoplot(xlim = c(20, 70), ylim = c(15, 34), dlat = 10, dlon = 10)
#' geoworld(database = "world", fill = TRUE, col = 30) #
#' 
#' geoplot(xlim = c(-10, 70), ylim = c(71, 81), b0 = 80, 
#'   dlat = 2, dlon = 10) # 0 must be high here else
#' geoworld(database = "world", fill = TRUE, col = 30) #the plot fails.  
#' 
#' # Lambert projection, 
#' geoplot(xlim = c(-10, 70), ylim = c(71, 81), 
#'   dlat = 2, dlon = 10, projection = "Lambert")
#' geoworld(database = "world", fill = TRUE, col = 30) 
#' }
#' @export geoworld
geoworld <-
function(regions = ".", exact = FALSE, boundary = TRUE, fill = FALSE, color = 1, lwd = 1, 
	lty = 1, plot = TRUE, type = "l", pch = ".", database = "world", 
	return.data = FALSE, outside=FALSE, allowed.size = 80000)
{
  geopar <- getOption("geopar")
  resolution <- 0 #1
  interior <- FALSE
  r <- 1.2
  doproj <- FALSE
 if(fill)
    interior <- TRUE
  if(geopar$projection == "Lambert") {
                                        # complicated borders in lat, lon
    p1 <- list(x = c(geopar$limx[1], mean(geopar$limx), geopar$limx[
                 1], geopar$limx[2]), y = c(geopar$limy[1], geopar$limy[
                                        2], geopar$limy[2], geopar$limy[2]))
    limits <- invProj(p1$x, p1$y, geopar$scale, geopar$b0, geopar$
                      b1, geopar$l1, geopar$projection)
    xlim <- c(limits$lon[3], limits$lon[4])
    ylim <- c(limits$lat[1], limits$lat[2])
  }
  else {
    limits <- invProj(geopar$limx, geopar$limy, geopar$scale, 
                      geopar$b0, geopar$b1, geopar$l1, geopar$projection)
    xlim <- c(limits$lon[1], limits$lon[2])
    ylim <- c(limits$lat[1], limits$lat[2])
  }
  xlim <- mean(xlim) + r * (xlim - mean(xlim))	# add
  ylim <- mean(ylim) + r * (ylim - mean(ylim))	# to get everything
                                        # parameter checks
                                        # turn the region names into a list of polygon numbers
  gon <- maps:::mapname(database, regions, exact)
  n <- length(gon)
  if(n == 0) stop("nothing to draw: no recognized region names")	
                                        # turn the polygon numbers into a list of polyline numbers
  line <- maps:::mapgetg(database, gon, fill, c(-1000000, 1000000), c(-1000000, 
                                                               1000000))
  if(length(line$number) == 0) stop(
             "nothing to draw: all regions out of bounds")	
                                        # turn the polyline numbers into x and y coordinates
  if(fill)
    coord <- maps:::mapgetl(database, line$number, c(-1000000, 1000000), c(
                                                                    -1000000, 1000000))
  else {
    l <- abs(line$number)
    if(boundary && interior)
      l <- unique(l)
    else if(boundary)
      l <- l[!match(l, l[duplicated(l)], FALSE)]
    else l <- l[duplicated(l)]
    coord <- maps:::mapgetl(database, l, xlim, ylim)
    if(length(coord) == 0)
      stop("all data out of bounds")
  }
  if(doproj) {
    coord <- Proj(coord$y, coord$x, geopar$scale, geopar$b0, geopar$
                  b1, geopar$l1, geopar$projection)
    coord$error <- FALSE
  }
                                        # for filled regions, turn NA breaks at polylines into
                                        # NA breaks at polygons, deleting polygons for which
                                        # there is a corresponding NA color
  if(fill) {
    gonsize <- line$size
    color <- rep(color, length = length(gonsize))
    keep <- !is.na(color)
    coord[c("x", "y")] <- maps:::makepoly(coord, gonsize, keep)
    color <- color[keep]
  }
  if(return.data) return(data.frame(lat = coord$y, lon = coord$x))	
                                        # do the plotting, if requested
  data <- data.frame(lat = coord$y, lon = coord$x)
  if(plot){
    if(fill) geopolygon(data, col = color, allowed.size = allowed.size)
    if(!fill) geolines(data, col = color, lwd = lwd, lty = lty)
  }
  return(invisible())
}

