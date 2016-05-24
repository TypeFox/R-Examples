# Envelope information from R spatial objects

#' Envelope of Spatial Objects
#'
#' Takes a \code{Spatial} object and returns a data frame with six columns
#' representing the envelope of each element in the \code{Spatial} object.
#'
#' @param obj object inheriting class \code{\link[sp:Spatial-class]{Spatial}}.
#' @details \code{obj} may be an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}},
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}},
#'   \code{\link[sp:SpatialLines-class]{SpatialLines}},
#'   \code{\link[sp:SpatialLinesDataFrame-class]{SpatialLinesDataFrame}},
#'   \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}, or
#'   \code{\link[sp:SpatialPolygonsDataFrame-class]{SpatialPolygonsDataFrame}},
#'   or a \code{list} in which each element is an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}} or
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.
#' @param centerfun function to apply to the x-axis limits and y-axis limits
#'   of the bounding box to obtain the x-coordinate and y-coordinate of the
#'   center of the bounding box.
#' @return A data frame with six columns named XMax, XMin, YMax, YMin, XCenter,
#'   and YCenter. The first four columns represent the corners of the bounding
#'   box of each element in \code{obj}. The last two columns represent the
#'   center of the bounding box of each element in \code{obj}. The number of
#'   rows in the returned data frame is the same as the length of the argument
#'   \code{obj}.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the columns of the returned data frame have the
#'   SpotfireColumnMetaData attribute set to enable TIBCO Spotfire to recognize
#'   them as containing envelope information.
#' @seealso Example usage at \code{\link{writeWKB}}
#' @export
writeEnvelope <- function(obj, centerfun = mean) {
  if(inherits(obj, c("SpatialPoints", "SpatialPointsDataFrame"), which = FALSE)) {

    SpatialPointsEnvelope(obj)

  } else if(inherits(obj, "list") && length(obj) > 0 &&
            all(vapply(
              X = obj,
              FUN = inherits,
              FUN.VALUE = logical(1),
              c("SpatialPoints", "SpatialPointsDataFrame"))
            )
  ) {

    ListOfSpatialPointsEnvelope(obj, centerfun = mean)

  } else if(inherits(obj, c("SpatialLines", "SpatialLinesDataFrame"), which = FALSE)) {

    SpatialLinesEnvelope(obj, centerfun = mean)

  } else if(inherits(obj, c("SpatialPolygons", "SpatialPolygonsDataFrame"), which = FALSE)) {

    SpatialPolygonsEnvelope(obj)

  } else {

    stop("obj must be an object of class SpatialPoints, SpatialPointsDataFrame, ",
         "SpatialLines, SpatialLinesDataFrame, SpatialPolygons, ",
         "or SpatialPolygonsDataFrame, or a list of objects of class ",
         "SpatialPoints or SpatialPointsDataFrame")

  }
}
