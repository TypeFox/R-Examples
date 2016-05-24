# Convert a SpatialPoints or SpatialPointsDataFrame object
#  to a well-known binary (WKB) geometry representation of points

#' Convert List of SpatialPoints to \acronym{WKB} MultiPoint
#'
#' Converts a list of objects of class \code{SpatialPoints} or
#' \code{SpatialPointsDataFrame} to a list of well-known binary (\acronym{WKB})
#' geometry representations of type MultiPoint.
#'
#' This function is called by the \code{\link{writeWKB}} function. Call the
#' \code{\link{writeWKB}} function instead of calling this function directly.
#'
#' Use this function when each item in the \acronym{WKB} representation should
#' represent multiple points. Use \code{\link{SpatialPointsToWKBPoint}} when
#' each item in the \acronym{WKB} representation should represent only one
#' point.
#'
#' @param obj a \code{list} in which each element is an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}} or
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.
#' @param endian The byte order (\code{"big"} or \code{"little"}) for encoding
#'   numeric types. The default is \code{"little"}.
#' @return A \code{list} with class \code{AsIs}. The length of the returned list
#'   is the same as the length of the argument \code{obj}. Each element of the
#'   returned list is a \code{\link[base]{raw}} vector consisting of a
#'   well-known binary (\acronym{WKB}) geometry representation of type
#'   MultiPoint.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the return value has the SpotfireColumnMetaData attribute
#'   set to enable TIBCO Spotfire to recognize it as a \acronym{WKB} geometry
#'   representation.
#' @examples
#' # load package sp
#' library(sp)
#'
#' # create a list of objects of class SpatialPoints
#' x1 = c(1, 2, 3, 4, 5)
#' y1 = c(3, 2, 5, 1, 4)
#' x2 <- c(9, 10, 11, 12, 13)
#' y2 <- c(-1, -2, -3, -4, -5)
#' Sp1 <- SpatialPoints(data.frame(x1, y1))
#' Sp2 <- SpatialPoints(data.frame(x2, y2))
#' obj <- list("a"=Sp1, "b"=Sp2)
#'
#' # convert to WKB MultiPoint
#' wkb <- wkb:::ListOfSpatialPointsToWKBMultiPoint(obj)
#'
#' # use as a column in a data frame
#' ds <- data.frame(ID = names(obj), Geometry = wkb)
#'
#' # calculate envelope columns and cbind to the data frame
#' coords <- wkb:::ListOfSpatialPointsEnvelope(obj)
#' ds <- cbind(ds, coords)
#' @seealso \code{\link{writeWKB}}, \code{\link{SpatialPointsToWKBPoint}},
#'   \code{\link{ListOfSpatialPointsEnvelope}}
#' @noRd
ListOfSpatialPointsToWKBMultiPoint <- function(obj, endian) {
  wkb <- lapply(X = obj, FUN = function(mypoints) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(4L, rc, size = 4, endian = endian)
    coords <- mypoints@coords
    writeBin(nrow(coords), rc, size = 4, endian = endian)
    apply(X = coords, MARGIN = 1, FUN = function(coord) {
      if(endian == "big") {
        writeBin(as.raw(0L), rc)
      } else {
        writeBin(as.raw(1L), rc)
      }
      writeBin(1L, rc, size = 4, endian = endian)
      writeBin(coord[1], rc, size = 8, endian = endian)
      writeBin(coord[2], rc, size = 8, endian = endian)
      NULL
    })
    rawConnectionValue(rc)
  })
  if(identical(version$language, "TERR")) {
    attr(wkb, "SpotfireColumnMetaData") <-
      list(ContentType = "application/x-wkb", MapChart.ColumnTypeId = "Geometry")
  }
  I(wkb)
}

#' Convert SpatialPoints to \acronym{WKB} Point
#'
#' Converts an object of class \code{SpatialPoints} or
#' \code{SpatialPointsDataFrame} to a list of well-known binary (\acronym{WKB})
#' geometry representations of type Point.
#'
#' This function is called by the \code{\link{writeWKB}} function. Call the
#' \code{\link{writeWKB}} function instead of calling this function directly.
#'
#' Use this function when each item in the \acronym{WKB} representation should
#' represent only one point. Use
#' \code{\link{ListOfSpatialPointsToWKBMultiPoint}} when each item in the
#' \acronym{WKB} representation should represent multiple points.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}} or
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.
#' @param endian The byte order (\code{"big"} or \code{"little"}) for encoding
#'   numeric types. The default is \code{"little"}.
#' @return A \code{list} with class \code{AsIs}. The length of the returned list
#'   is the same as the length of the argument \code{obj}. Each element of the
#'   returned list is a \code{\link[base]{raw}} vector consisting of a
#'   well-known binary (\acronym{WKB}) geometry representation of type Point.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the return value has the SpotfireColumnMetaData attribute
#'   set to enable TIBCO Spotfire to recognize it as a \acronym{WKB} geometry
#'   representation.
#' @examples
#' # create an object of class SpatialPoints
#' x = c(1, 2, 3, 4, 5)
#' y = c(3, 2, 5, 1, 4)
#' Sp <- SpatialPoints(data.frame(x, y))
#'
#' # convert to WKB Point
#' wkb <- wkb:::SpatialPointsToWKBPoint(Sp)
#'
#' # use as a column in a data frame
#' ds <- data.frame(ID = c("a", "b", "c", "d", "e"), Geometry = wkb)
#'
#' # calculate envelope and center columns and cbind to the data frame
#' coords <- wkb:::SpatialPointsEnvelope(Sp)
#' ds <- cbind(ds, coords)
#' @seealso \code{\link{writeWKB}},
#'   \code{\link{ListOfSpatialPointsToWKBMultiPoint}},
#'   \code{\link{SpatialPointsEnvelope}}
#' @noRd
SpatialPointsToWKBPoint <- function(obj, endian) {
  wkb <- lapply(apply(X = obj@coords, MARGIN = 1, FUN = function(coord) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(1L, rc, size = 4, endian = endian)
    writeBin(coord[1], rc, size = 8, endian = endian)
    writeBin(coord[2], rc, size = 8, endian = endian)
    list(rawConnectionValue(rc))
  }), unlist)
  if(identical(version$language, "TERR")) {
    attr(wkb, "SpotfireColumnMetaData") <-
      list(ContentType = "application/x-wkb", MapChart.ColumnTypeId = "Geometry")
  }
  I(wkb)
}

#' Envelope of List of SpatialPoints
#'
#' Takes a list of objects of class \code{SpatialPoints} or
#' \code{SpatialPointsDataFrame} and returns a data frame with six columns
#' representing the envelope of each object of class \code{SpatialPoints} or
#' \code{SpatialPointsDataFrame}.
#'
#' This function is called by the \code{\link{writeEnvelope}} function. Call the
#' \code{\link{writeEnvelope}} function instead of calling this function
#' directly.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}} or
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.
#' @param centerfun function to apply to the x-axis limits and y-axis limits of
#'   the bounding box to obtain the x-coordinate and y-coordinate of the center
#'   of the bounding box.
#' @return A data frame with six columns named XMax, XMin, YMax, YMin, XCenter,
#'   and YCenter. The first four columns represent the corners of the bounding
#'   box of each object of class \code{SpatialPoints} or
#'   \code{SpatialPointsDataFrame}. The last two columns represent the center of
#'   the bounding box of each object of class \code{SpatialPoints} or
#'   \code{SpatialPointsDataFrame}. The number of rows in the returned data
#'   frame is the same as the length of the argument \code{obj}.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the columns of the returned data frame have the
#'   SpotfireColumnMetaData attribute set to enable TIBCO Spotfire to recognize
#'   them as containing envelope information.
#' @seealso \code{\link{writeEnvelope}}
#'
#'   Example usage at \code{\link{ListOfSpatialPointsToWKBMultiPoint}}
#' @noRd
#' @importFrom sp bbox
ListOfSpatialPointsEnvelope <- function(obj, centerfun = mean) {
  if(is.character(centerfun)) {
    centerfun <- eval(parse(text = centerfun))
  }
  coords <- as.data.frame(t(vapply(X = obj, FUN = function(mypoints) {
    c(XMax = bbox(mypoints)[1, "max"],
      XMin = bbox(mypoints)[1, "min"],
      YMax = bbox(mypoints)[2, "max"],
      YMin = bbox(mypoints)[2, "min"],
      XCenter = centerfun(bbox(mypoints)[1, ], na.rm = TRUE),
      YCenter = centerfun(bbox(mypoints)[2, ], na.rm = TRUE))
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

#' Envelope of SpatialPoints
#'
#' Takes an object of class \code{SpatialPoints} or
#' \code{SpatialPointsDataFrame} and returns a data frame with six columns
#' representing the envelope of each point (which is each point itself).
#'
#' This function is called by the \code{\link{writeEnvelope}} function. Call the
#' \code{\link{writeEnvelope}} function instead of calling this function
#' directly.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}} or
#'   \code{\link[sp:SpatialPointsDataFrame-class]{SpatialPointsDataFrame}}.
#' @return A data frame with six columns named XMax, XMin, YMax, YMin, XCenter,
#'   and YCenter. The first four columns represent the corners of the bounding
#'   box of each point. The last two columns represent the center of the
#'   bounding box of each point. (Note that the bounding box of a point is the
#'   point itself.) The number of rows in the returned data frame is the same as
#'   the length of the argument \code{obj}.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the columns of the returned data frame have the
#'   SpotfireColumnMetaData attribute set to enable TIBCO Spotfire to recognize
#'   them as containing envelope information.
#' @seealso \code{\link{writeEnvelope}}
#'
#'   Example usage at \code{\link{SpatialPointsToWKBPoint}}
#' @noRd
SpatialPointsEnvelope <- function(obj) {
  coords <- as.data.frame(t(apply(X = obj@coords, MARGIN = 1, FUN = function(coord) {
    c(coord[1], # XMax
      coord[1], # XMin
      coord[2], # YMax
      coord[2], # YMin
      coord[1], # XCenter
      coord[2]) # YCenter
  })))
  colnames(coords) <- c("XMax", "XMin", "YMax", "YMin", "XCenter", "YCenter")
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
