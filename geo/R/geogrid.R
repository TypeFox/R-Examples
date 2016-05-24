#' Plots a grid.
#' 
#' Plots a grid defined by the vectors lat, lon. The grid is plotted on a graph
#' initialized by geoplot.  lon gives the meridians plotted and lat the
#' parallels plotted.
#' 
#' 
#' @param lat,lon Latitude and longitude of data ( or x and y coordinates),
#' negative for southern latitudes and western longitudes.  May be supplied as
#' two vectors or as a dataframe lat (or x) including vectors \code{lat$lat}
#' and \code{lat$lon} (\code{x$x} and \code{x$y} if projection = none).
#' @param col Color number used, default value is 1 (black).
#' @param type "l" means line and "p" points.  Default is "l".
#' @param lwd Linewidth.  Default value is the value set when the program was
#' called.
#' @param lty Linetype.  Default value is the value set when the program was
#' called.
#' @param pch Type of symbol at gridpoints default is "+".
#' @param nx sets smoothness of curved Lambert parallels
#' @return No values returned.
#' @seealso \code{\link{geoplot}}, \code{\link{geolines}},
#' \code{\link{geopolygon}}, \code{\link{geotext}}, \code{\link{geosymbols}},
#' \code{\link{geopar}}, \code{\link{geolocator}}, \code{\link{geocontour}}.
#' @examples
#' 
#' \dontrun{       geogrid(latgr, longgr)
#' 
#'        codgrd <- list(lat = seq(62, 68, by = 0.1), lon = seq(-28, -10, 0.25))
#'        geogrid(codgrd)   # a fine grid of Iceland and neighbouring seas
#'        geoplot(new = T)   
#' }
#' @export geogrid
geogrid <-
function (lat, lon = 0, col = 1, type = "l", lwd = 0, lty = 0, 
    pch = "+", nx = 5) 
{
    geopar <- getOption("geopar")
    oldpar <- selectedpar()
    if (length(lon) == 1) {
        if (geopar$projection == "none") {
            lon <- lat$y
            lat <- lat$x
        }
        else {
            lon <- lat$lon
            lat <- lat$lat
        }
    }
    if (geopar$projection == "Lambert") 
        nx <- nx
    else nx <- 1
    if (geopar$projection != "none") {
        if (mean(lat, na.rm = T) > 1000) {
            lat <- geoconvert(lat)
            lon <- -geoconvert(lon)
        }
    }
    if (type == "l") {
        llon <- length(lon)
        llat <- length(lat)
        latgr <- t(matrix(lat, llat, llon))
        longr <- matrix(lon, llon, llat)
        latgr <- rbind(latgr, rep(NA, ncol(latgr)))
        longr <- rbind(longr, rep(NA, ncol(longr)))
        geolines(latgr, longr, 
          col = col, lwd = lwd, lty = lty, nx = nx)
        llon <- length(lon)
        llat <- length(lat)
        latgr <- matrix(lat, llat, llon)
        longr <- t(matrix(lon, llon, llat))
        latgr <- rbind(latgr, rep(NA, ncol(latgr)))
        longr <- rbind(longr, rep(NA, ncol(longr)))
        geolines(latgr, longr, 
          col = col, lwd = lwd, lty = lty, nx = nx)
    }
    else {
        llon <- length(lon)
        llat <- length(lat)
        latgr <- c(t(matrix(lat, llat, llon)))
        longr <- c(matrix(lon, llon, llat))
        geopoints(latgr, longr, pch = pch)
    }
    return(invisible())
    par(oldpar)
}
