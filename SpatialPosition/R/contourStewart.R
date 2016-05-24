#' @title Create a SpatialPolygonsDataFrame or a SpatialLinesDataFrame from a 
#' Stewart Raster
#' @name contourStewart
#' @description 
#' \code{contourStewart} is deprecated. \cr
#' To obtain contour lines use \code{\link{rasterToContour}} from raster package. \cr
#' To obtain contour polygons use \code{\link{rasterToContourPoly}} from SpatialPosition package. \cr\cr\cr
#' This function creates a SpatialPolygonsDataFrame or SpatialLinesDataFrame contour from the Stewart raster.
#' @param x raster; output of the \code{\link{rasterStewart}} function. The raster must contain only positive values.
#' @param breaks numeric; a vector of break values. 
#' @param mask SpatialPolygonsDataFrame; mask used to clip contour shapes.
#' @param type character; "poly" or "line". WARNING: the poly option is experimental (see details). It needs the rgeos package.
#' @return The ouput of the function is a SpatialPolygonsDataFrame (\code{type = "poly"}) 
#' or a SpatialLinesDataFrame (\code{type = "line"}). 
#' The data frame of the outputed SpatialPolygonsDataFrame contains four fields: 
#' id (id of each polygon), min and max (minimum and maximum breaks of the polygon), 
#' mean (center value of the class)
#' @details To obtain a correct SpatialPolygonsDataFrame of potentials follow theses steps: \itemize{
#' \item{Step 1: Create a SpatialPointsDataFrame of potentials with the 
#' stewart function. Do not enter an unknownpts layer, set a resolution, 
#' and set a SpatialPolygonsDataFrame (spmask) as mask.}
#' \item{Step 2: Create a raster from the SpatialPointsDataFrame of potentials 
#' with the rasterStewart function without using a mask.}
#' \item{Step 3: Create the SpatialPolygonsDataFrame of potentials with the 
#' contourStewart function and use the same spmask SpatialPolygonsDataFrame (Step1) as mask.}
#' }
#' See also the second example in the examples section.
#' @seealso \link{stewart}, \link{rasterStewart}, \link{plotStewart}, 
#' \link{quickStewart}, \link{CreateGrid}, \link{CreateDistMatrix}.
#' @import sp
#' @import raster
#' @examples
#' data("spatData")
#' \dontrun{
#' #### Example with type = "line"
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE,
#'                      mask = spatMask)
#' # Create a raster of potentials values
#' mystewartraster <- rasterStewart(x = mystewart, mask = spatMask)
#' # Display the raster and get break values
#' break.values <- plotStewart(x = mystewartraster)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourpoly <- contourStewart(x = mystewartraster,
#'                                        breaks = break.values,
#'                                        type = "line")
#' # Display the Map
#' plot(spatMask, add=TRUE)
#' plot(mystewartcontourpoly, border = "grey40",add = TRUE)
#' plot(spatPts, cex = 0.8, pch = 20, col  = "black", add = TRUE)
#' 
#' 
#' 
#' #### Example with type = "poly"
#' mystewart <- stewart(knownpts = spatPts, varname = "Capacite",
#'                      typefct = "exponential", span = 1000, beta = 3,
#'                      resolution = 50, longlat = FALSE,
#'                      mask = spatMask)
#' # Create a raster of potentials valuesn, no mask
#' mystewartraster <- rasterStewart(x = mystewart)
#' # Display the raster and get break values
#' break.values <- plotStewart(x = mystewartraster)
#' # Create contour SpatialLinesDataFrame
#' mystewartcontourpoly <- contourStewart(x = mystewartraster,
#'                                        breaks = break.values,
#'                                        mask = spatMask,
#'                                        type = "poly")
#' # Display the map
#' library(cartography)
#' opar <- par(mar = c(0,0,1.1,0))
#' choroLayer(spdf = mystewartcontourpoly, 
#'            df = mystewartcontourpoly@data, 
#'            var = "mean", legend.pos = "topleft",
#'            breaks = break.values, border = "grey90", 
#'            lwd = 0.2, 
#'            legend.title.txt = "Potential number\nof beds in the\nneighbourhood", 
#'            legend.values.rnd = 0)
#' plot(spatMask, add = TRUE)
#' propSymbolsLayer(spdf = spatPts, df = spatPts@data, var = "Capacite",
#'                  legend.title.txt = "Number of beds", 
#'                  col = "#ff000020")
#' layoutLayer(title = "Global Accessibility to Public Hospitals", 
#'             south = TRUE, sources = "", author = "")
#' par(opar)
#' }
#' @export
contourStewart <- function(x, breaks, mask, type = "line"){
  if (type=="line"){
    .Deprecated(new = "rasterToContour", package = "raster") 
    return(rasterToContour(x = x, levels = breaks[1:(length(breaks)-1)]))
  } 
  if (type=="poly"){
    .Deprecated("rasterToContourPoly", package = "SpatialPosition")
    return(rasterToContourPoly(r = x, breaks = breaks, mask = mask))
  }
}