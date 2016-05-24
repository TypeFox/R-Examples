#' @title Plot Samples on World Map
#' @description Plot a set of samples on a world map.
#' 
#' @param lon,lat vectors giving the longitude and latitude of points to plot.
#' @param lon.range,lat.range vectors giving the minimum and maximum longitude 
#'   and latitude of the map. If the first value in \code{lon.range} is greater 
#'   than the second value, then a Pacific-centric map (\code{\link[mapdata]{world2Hires}}) 
#'   is used and continents will not be filled in.
#' @param main main title for the plot.
#' @param pch point symbol to use.
#' @param pt.cex point size to use.
#' @param col point color.
#' @param bg background color of point.
#' @param n,lon.n,lat.n the number of tick marks desired. Can be specified 
#'   separately for longitude (\code{lon.n}) or latitude (\code{lat.n}). 
#'   See \code{\link{pretty}} for more details.
#' 
#' @return original \code{\link{par}} settings for \code{mar} and \code{oma}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{lat.lon.axes}}
#' 
#' @examples
#' 
#' # Some random points around San Diego, CA
#' lat <- runif(30, 32.5, 33.3)
#' lon <- runif(30, -118.3, -117.5)
#' lat.range <- c(32.4, 33.6)
#' lon.range <- c(-118.6, -117)
#' sample.map(lat, lon, lat.range, lon.range)
#' 
#' # Some random points in the Pacific
#' lat <- runif(30, 20, 50)
#' lon <- c(runif(15, 150, 180), runif(15, -180, -120))
#' lat.range <- c(10, 75)
#' lon.range <- c(110, -110)
#' sample.map(lat, lon, lat.range, lon.range)
#' 
#' @importFrom stats sd
#' @importFrom graphics par points mtext box
#' @importFrom maps map
#' @import mapdata
#' @export sample.map
#' 
sample.map <- function(lat, lon, lat.range = NULL, lon.range = NULL, main = NULL, 
                       pch = 19, pt.cex = 1, col = "black", bg = col, n = 5, 
                       lon.n = n, lat.n = n) {
  if(is.null(lon.range)) lon.range <- range(lon)
  if(is.null(lat.range)) lat.range <- range(lat)
  pacific.cent <- lon.range[1] > lon.range[2]
  if(pacific.cent) {
    lon <- ifelse(lon < 0, 360 + lon, lon)
    lon.range <- ifelse(lon.range < 0, 360 + lon.range, lon.range)
  }
  
  has.loc <- !is.na(lon) & !is.na(lat) 
  in.lon.range <- lon >= min(lon.range) & lon <= max(lon.range)
  in.lat.range <- lat >= min(lat.range) & lat <= max(lat.range)
  to.plot <- has.loc & in.lon.range & in.lat.range

  if(!is.null(main)) main <- paste(main, " (n = ", sum(to.plot), ")", sep = "")
  if(length(pch) == length(lon)) pch <- pch[to.plot]
  if(length(pt.cex) == length(lon)) pt.cex <- pt.cex[to.plot]
  if(length(col) == length(lon)) col <- col[to.plot]

  op <- par(mar = c(3, 5, ifelse(is.null(main), 3, 5), 5) + 0.1, oma = c(1, 1, 1, 1))
    
  if(pacific.cent) {
    map("world2Hires", xlim = lon.range, ylim = lat.range)
  } else {
    map("worldHires", fill = TRUE, col = "wheat3", xlim = lon.range, ylim = lat.range)
  }
  
  points(lon[to.plot], lat[to.plot], pch = pch, cex = pt.cex, col = col, bg = bg)

  lat.lon.axes(n = n)
  
  if(!is.null(main)) mtext(main, line = 3, cex = 1.5)
  box(lwd = 2)
  
  invisible(op)
}