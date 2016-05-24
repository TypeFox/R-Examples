#' Create an anomaly RasterStack 
#' 
#' @description
#' The function creates an anomaly RasterStack either based on the
#' overall mean of the original stack, or a supplied reference RasterLayer.
#' For the creation of seasonal anomalies use \code{\link{deseason}}.
#' 
#' @param x a RasterStack
#' @param reference an optional RasterLayer to be used as the reference 
#' @param ... additional arguments passed to \code{\link{calc}} 
#' which is used under the hood
#' 
#' @return an anomaly RasterStack
#' 
#' @seealso
#' \code{\link{deseason}}, \code{\link{denoise}}, \code{\link{calc}}
#' 
#' @export anomalize
#' 
#' @examples
#' data(australiaGPCP)
#' 
#' aus_anom <- anomalize(australiaGPCP)
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(australiaGPCP[[10]], main = "original")
#' plot(aus_anom[[10]], main = "anomalized")
#' par(opar)
anomalize <- function(x, 
                      reference = NULL, 
                      ...) {
  
  if (is.null(reference)) {
    mn <- raster::calc(x, fun = mean, ...)
  } else {
    mn <- reference
  }
  
  return(x - mn)
  
}