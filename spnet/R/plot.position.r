#' Plot a map labelled with the ID numbering
#'
#' The \code{graph.map.plot.position} function allows to plot maps defined as for example \code{SpatialNetwork} or \code{SpatialPolygons} objects, and render the ID numbering.
#'
#' @param x an object for which a \code{graph.map.plot.position} method is defined.
#' @param label a character of length 1 for prefixing seat numbering.
#' @param ... other arguments to pass to the plot function. The main usage is setting the \code{cex} value.
#'
#' @return NULL
#' 
#' @family plot map
#' @importFrom graphics text
#' @export
#' @examples
#' ## The world map
#' data(world.map.simplified, package = "spnet")
#'
#' graph.map.plot.position(world.map.simplified)
#' graph.map.plot.position(world.map.simplified, cex = 0.4)
#' graph.map.plot.position(world.map.simplified, label = 'ID ', cex = 0.3)
setGeneric(
  'graph.map.plot.position',
  function(
    x,
    label = '',
    ...
  ) {
    standardGeneric("graph.map.plot.position")
  }
)

#' @describeIn graph.map.plot.position method for \code{SpatialPolygons} objects.
setMethod(
  f = 'graph.map.plot.position',
  signature = 'SpatialPolygons',
  definition = function(x, label, ...) {
    if(length(x) == 0)
      stop("The map is empty. Please define a valid map.")
    plot(x, ...)
    text(coordinates(x), labels=paste(label, row.names(coordinates(x)), sep=''), ...)
  }
)

#' @describeIn graph.map.plot.position method for \code{SpatialNetwork} objects.
setMethod(
  f = 'graph.map.plot.position',
  signature = 'SpatialNetwork',
  definition = function(x, label, ...) {
    if(length(x@map) == 0)
      stop("The map is empty. Please define a valid map.")
    getMethod('graph.map.plot.position', 'SpatialPolygons')(x@map, label, ...)
  }
)