#' @title Latitude and Longitude axes
#' @description Add latitude and longitude axes to a map.
#' 
#' @param n,lon.n,lat.n the number of tick marks desired. Can be specified separately for longitude (\code{lon.n}) 
#' or latitude (\code{lat.n}). See \code{\link{pretty}} for more details.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom graphics par axis
#' @export
#' 
lat.lon.axes <- function(n = 5, lon.n = n, lat.n = n) {
  lon.range <- par("usr")[1:2]
  lon <- list(ticks = pretty(lon.range, n = lon.n))
  lon$labels <- parse(text = sapply(lon$ticks, function(i) {
    i <- ifelse(i > 180, i - 360, i)
    i <- ifelse(i < -180, 360 - i, i)
    a <- ifelse(i < 0, -1 * i, i)
    b <- ifelse(i < 0, "~W", "~E")
    b <- ifelse(a %in% c(0, 180), "", b)
    paste(a, "*degree", b, sep = "")
  }))
  
  lat.range <- par("usr")[3:4]
  lat <- list(ticks = pretty(lat.range, n = lat.n))
  lat$labels <- parse(text = sapply(lat$ticks, function(i) {
    a <- ifelse(i < 0, -1 * i, i)
    b <- ifelse(i < 0, "~S", "~N")
    b <- ifelse(a %in% c(0, 90), "", b)
    paste(a, "*degree", b, sep = "")
  }))
  
  axis(1, at = lon$ticks, labels = lon$labels)
  axis(2, at = lat$ticks, labels = lat$labels, las = 1)
  axis(3, at = lon$ticks, labels = lon$labels)
  axis(4, at = lat$ticks, labels = lat$labels, las = 1)
  invisible(NULL)
}
