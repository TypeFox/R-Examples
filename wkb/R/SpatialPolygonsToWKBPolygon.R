# Convert a SpatialPolygons or SpatialPolygonsDataFrame object
#  to a well-known binary (WKB) geometry representation of polygons

#' Convert SpatialPolygons to \acronym{WKB} Polygon
#'
#' Converts an object of class \code{SpatialPolygons} or
#' \code{SpatialPolygonsDataFrame} to a list of well-known binary
#' (\acronym{WKB}) geometry representations of type Polygon.
#'
#' This function is called by the \code{\link{writeWKB}} function. Call the
#' \code{\link{writeWKB}} function instead of calling this function directly.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}} or
#'   \code{\link[sp:SpatialPolygonsDataFrame-class]{SpatialPolygonsDataFrame}}.
#' @param endian The byte order (\code{"big"} or \code{"little"}) for encoding
#'   numeric types. The default is \code{"little"}.
#' @return A \code{list} with class \code{AsIs}. The length of the returned list
#'   is the same as the length of the argument \code{obj}. Each element of the
#'   returned list is a \code{\link[base]{raw}} vector consisting of a
#'   well-known binary (\acronym{WKB}) geometry representation of type Polygon.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the return value has the SpotfireColumnMetaData attribute
#'   set to enable TIBCO Spotfire to recognize it as a \acronym{WKB} geometry
#'   representation.
#' @examples
#' # load package sp
#' library(sp)
#'
#' # create an object of class SpatialPolygons
#' triangle <- Polygons(
#'  list(
#'    Polygon(data.frame(x = c(2, 2.5, 3, 2), y = c(2, 3, 2, 2)))
#'  ), "triangle")
#' rectangles <- Polygons(
#'   list(
#'     Polygon(data.frame(x = c(0, 0, 1, 1, 0), y = c(0, 1, 1, 0, 0))),
#'     Polygon(data.frame(x = c(0, 0, 2, 2, 0), y = c(-2, -1, -1, -2, -2)))
#'   ), "rectangles")
#' Sp <- SpatialPolygons(list(triangle, rectangles))
#'
#' # convert to WKB Polygon
#' wkb <- wkb:::SpatialPolygonsToWKBPolygon(Sp)
#'
#' # use as a column in a data frame
#' ds <- data.frame(ID = names(Sp), Geometry = wkb)
#'
#' # calculate envelope columns and cbind to the data frame
#' coords <- wkb:::SpatialPolygonsEnvelope(Sp)
#' ds <- cbind(ds, coords)
#' @seealso \code{\link{writeWKB}}, \code{\link{SpatialPolygonsEnvelope}}
#' @noRd
SpatialPolygonsToWKBPolygon <- function(obj, endian) {
  wkb <- lapply(X = obj@polygons, FUN = function(mypolygon) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(3L, rc, size = 4, endian = endian)
    rings <- mypolygon@Polygons
    writeBin(length(rings), rc, size = 4, endian = endian)
    lapply(X = rings, FUN = function(ring) {
      coords <- ring@coords
      writeBin(nrow(coords), rc, size = 4, endian = endian)
      apply(X = coords, MARGIN = 1, FUN = function(coord) {
        writeBin(coord[1], rc, size = 8, endian = endian)
        writeBin(coord[2], rc, size = 8, endian = endian)
        NULL
      })
    })
    rawConnectionValue(rc)
  })
  if(identical(version$language, "TERR")) {
    attr(wkb, "SpotfireColumnMetaData") <-
      list(ContentType = "application/x-wkb", MapChart.ColumnTypeId = "Geometry")
  }
  I(wkb)
}

#' Envelope of SpatialPolygons
#'
#' Takes an object of class \code{SpatialPolygons} or
#' \code{SpatialPolygonsDataFrame} and returns a data frame with six columns
#' representing the envelope of each object of class \code{Polygons}.
#'
#' This function is called by the \code{\link{writeEnvelope}} function. Call the
#' \code{\link{writeEnvelope}} function instead of calling this function
#' directly.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}} or
#'   \code{\link[sp:SpatialPolygonsDataFrame-class]{SpatialPolygonsDataFrame}}.
#' @return A data frame with six columns named XMax, XMin, YMax, YMin, XCenter,
#'   and YCenter. The first four columns represent the corners of the bounding
#'   box of each object of class \code{Polygons}. The last two columns represent
#'   the center of the bounding box of each object of class \code{Polygons}. The
#'   number of rows in the returned data frame is the same as the length of the
#'   argument \code{obj}.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the columns of the returned data frame have the
#'   SpotfireColumnMetaData attribute set to enable TIBCO Spotfire to recognize
#'   them as containing envelope information.
#' @seealso \code{\link{writeEnvelope}}
#'
#'   Example usage at \code{\link{SpatialPolygonsToWKBPolygon}}
#' @noRd
#' @importFrom sp bbox
SpatialPolygonsEnvelope <- function(obj) {
  coords <- as.data.frame(t(vapply(X = obj@polygons, FUN = function(mypolygon) {
    c(XMax = bbox(mypolygon)["x", "max"],
      XMin = bbox(mypolygon)["x", "min"],
      YMax = bbox(mypolygon)["y", "max"],
      YMin = bbox(mypolygon)["y", "min"],
      XCenter = mypolygon@labpt[1],
      YCenter = mypolygon@labpt[2])
  }, FUN.VALUE = rep(0, 6))))
  if(identical(version$language, "TERR")) {
    attr(coords$XMax, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "XMax")
    attr(coords$XMin, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "XMin")
    attr(coords$YMax, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "YMax")
    attr(coords$YMin, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "YMin")
    attr(coords$XCenter, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "XCenter")
    attr(coords$YCenter, "SpotfireColumnMetaData") <- list(MapChart.ColumnTypeId = "YCenter")
  }
  coords
}
