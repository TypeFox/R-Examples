#' Deprecated functions in trip
#' 
#' These functions will be declared defunct in a future release.
#' 
#' @name trip.split.exact
#' @aliases trip-deprecated trip.split.exact as.trip.SpatialLinesDataFrame
#' tripTransform as.ltraj.trip
#' @param x see \code{\link{cut.trip}}
#' @param dates see \code{\link{cut.trip}}
#' @param from trip object
#' @param crs CRS object, or PROJ.4 string accepted by \code{\link[sp]{CRS}}
#' @param \dots Further arguments to \code{\link[rgdal]{spTransform}}
#' @seealso
#' 
#' \code{\link{cut.trip}}, \code{\link{as.Other}}
#' 
NULL

#' @rdname trip.split.exact
#' @export
as.SpatialLinesDataFrame.trip <- function (from) 
  {
    .Deprecated('as(x, "SpatialLinesDataFrame")')
    as(from, "SpatialLinesDataFrame")
  }

#' @rdname trip.split.exact
#' @export
trip.split.exact <- function(x, dates) {
  .Deprecated("cut.trip")
  cut(x, dates)
}

#' @rdname trip.split.exact
#' @param  xy \code{trip} object
#' @export
as.ltraj.trip <- function(xy) {
  .Deprecated('as(x, "ltraj")')
  as(xy, "ltraj")
}

##' @rdname trip.split.exact
##' @export
as.trip.SpatialLinesDataFrame <- function(from) {
  .Deprecated('as(x, "SpatialLinesDataFrame") or explode(x) ... the original definition was an error, there is no general coercion method available for SpatialLinesDataFrame=>trip')

  ##as.SpatialLinesDataFrame.trip(from)
  as(from, "SpatialLinesDataFrame")
}

## removed depends sp, suggests rgdal deprecate this, replace with
## spTransform method below

##' @rdname trip.split.exact
##' @export
tripTransform <- function(x, crs, ...) {
  .Deprecated("spTransform")
  if (! inherits(crs, "CRS")) crs <- CRS(crs)
  spTransform(x, crs, ...)
}
