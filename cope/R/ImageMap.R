#' Plots an image on a map of the world.
#'
#' @param lon Longitude.
#' @param lat Latitude.
#' @param img The image to plot as a matrix of dimensions 
#'            c(length(lon),length(lat)).
#' @param mask Matrix of dimensions c(length(lon),length(lat)) defining a mask 
#'             to cut out of the picture.
#' @param xlab Label for the x-axis passed to fields::image.plot. Default is "Longitude".
#' @param ylab Label for the y-axis passed to fields::image.plot. Default is "Latitude".
#' @param ... Additional graphical parameters passed to fields::image.plot.
#' @return NULL
ImageMap = function(lon, lat, img, mask=NULL, xlab='Longitude', ylab='Latitude', ...) {
  if(!is.null(mask)) {
    img = img*mask
    xlim = lon[range(which(rowSums(mask, na.rm=TRUE)>0))]
    ylim = lat[range(which(colSums(mask, na.rm=TRUE)>0))]
  } else {
    xlim = range(lon)
    ylim = range(lat)
  }
  fields::image.plot(lon, lat, img, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  maps::map('world', add=TRUE)
}