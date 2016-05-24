# Convert a SpatialLines or SpatialLinesDataFrame object
#  to a well-known binary (WKB) geometry representation of line segments

#' Convert SpatialLines to \acronym{WKB} MultiLineString
#'
#' Converts an object of class \code{SpatialLines} or
#' \code{SpatialLinesDataFrame} to a list of well-known binary (\acronym{WKB})
#' geometry representations of type MultiLineString.
#'
#' This function is called by the \code{\link{writeWKB}} function. Call the
#' \code{\link{writeWKB}} function instead of calling this function directly.
#'
#' The argument \code{obj} may have multiple objects of class \code{Lines} in
#' each position of the \code{list} in slot \code{lines}.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialLines-class]{SpatialLines}} or
#'   \code{\link[sp:SpatialLinesDataFrame-class]{SpatialLinesDataFrame}}.
#' @param endian The byte order (\code{"big"} or \code{"little"}) for encoding
#'   numeric types. The default is \code{"little"}.
#' @return A \code{list} with class \code{AsIs}. The length of the returned list
#'   is the same as the length of the argument \code{obj}. Each element of the
#'   returned list is a \code{\link[base]{raw}} vector consisting of a
#'   well-known binary (\acronym{WKB}) geometry representation of type
#'   MultiLineString.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the return value has the SpotfireColumnMetaData attribute
#'   set to enable TIBCO Spotfire to recognize it as a \acronym{WKB} geometry
#'   representation.
#' @examples
#' # load package sp
#' library(sp)
#'
#' # create an object of class SpatialLines
#' l1 <- data.frame(x = c(1, 2, 3), y = c(3, 2, 2))
#' l1a <- data.frame(x = l1[, 1] + .05, y = l1[, 2] + .05)
#' l2 <- data.frame(x = c(1, 2, 3), y = c(1, 1.5, 1))
#' Sl1 <- Line(l1)
#' Sl1a <- Line(l1a)
#' Sl2 <- Line(l2)
#' S1 <- Lines(list(Sl1, Sl1a), ID = "a")
#' S2 <- Lines(list(Sl2), ID = "b")
#' Sl <- SpatialLines(list(S1, S2))
#'
#' # convert to WKB MultiLineString
#' wkb <- wkb:::SpatialLinesToWKBMultiLineString(Sl)
#'
#' # use as a column in a data frame
#' ds <- data.frame(ID = names(Sl), Geometry = wkb)
#'
#' # calculate envelope columns and cbind to the data frame
#' coords <- wkb:::SpatialLinesEnvelope(Sl)
#' ds <- cbind(ds, coords)
#' @seealso \code{\link{writeWKB}}, \code{\link{SpatialLinesToWKBLineString}},
#'   \code{\link{SpatialLinesEnvelope}}
#' @noRd
SpatialLinesToWKBMultiLineString <- function(obj, endian) {
  wkb <- lapply(X = obj@lines, FUN = function(mylines) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(5L, rc, size = 4, endian = endian)
    lineStrings <- mylines@Lines
    writeBin(length(lineStrings), rc, size = 4, endian = endian)
    lapply(X = lineStrings, FUN = function(myline) {
      if(endian == "big") {
        writeBin(as.raw(0L), rc)
      } else {
        writeBin(as.raw(1L), rc)
      }
      writeBin(2L, rc, size = 4, endian = endian)
      coords <- myline@coords
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

#' Convert SpatialLines to \acronym{WKB} LineString
#'
#' Converts an object of class \code{SpatialLines} or
#' \code{SpatialLinesDataFrame} to a list of well-known binary (\acronym{WKB})
#' geometry representations of type LineString.
#'
#' The argument \code{obj} must have only one object of class \code{Lines} in
#' each position of the \code{list} in slot \code{lines}. If there are multiple
#' objects of class \code{Lines} in each position of the \code{list} in slot
#' \code{lines}, use \code{\link{SpatialLinesToWKBMultiLineString}}.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialLines-class]{SpatialLines}} or
#'   \code{\link[sp:SpatialLinesDataFrame-class]{SpatialLinesDataFrame}}.
#' @param endian The byte order (\code{"big"} or \code{"little"}) for encoding
#'   numeric types. The default is \code{"little"}.
#' @return A \code{list} with class \code{AsIs}. The length of the returned list
#'   is the same as the length of the argument \code{obj}. Each element of the
#'   returned list is a \code{\link[base]{raw}} vector consisting of a
#'   well-known binary (\acronym{WKB}) geometry representation of type
#'   LineString.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the return value has the SpotfireColumnMetaData attribute
#'   set to enable TIBCO Spotfire to recognize it as a \acronym{WKB} geometry
#'   representation.
#' @examples
#' # create an object of class SpatialLines
#' l1 <- data.frame(x = c(1, 2, 3), y = c(3, 2, 2))
#' l2 <- data.frame(x = c(1, 2, 3), y = c(1, 1.5, 1))
#' Sl1 <- Line(l1)
#' Sl2 <- Line(l2)
#' S1 <- Lines(list(Sl1), ID = "a")
#' S2 <- Lines(list(Sl2), ID = "b")
#' Sl <- SpatialLines(list(S1, S2))
#'
#' # convert to WKB LineString
#' wkb <- wkb:::SpatialLinesToWKBLineString(Sl)
#'
#' # use as a column in a data frame
#' ds <- data.frame(ID = names(Sl), Geometry = wkb)
#'
#' # calculate envelope columns and cbind to the data frame
#' coords <- wkb:::SpatialLinesEnvelope(Sl)
#' ds <- cbind(ds, coords)
#' @seealso \code{\link{writeWKB}},
#'   \code{\link{SpatialLinesToWKBMultiLineString}},
#'   \code{\link{SpatialLinesEnvelope}}
#' @noRd
SpatialLinesToWKBLineString <- function(obj, endian) {
  wkb <- lapply(X = obj@lines, FUN = function(mylines) {
    rc <- rawConnection(raw(0), "r+")
    on.exit(close(rc))
    if(endian == "big") {
      writeBin(as.raw(0L), rc)
    } else {
      writeBin(as.raw(1L), rc)
    }
    writeBin(2L, rc, size = 4, endian = endian)
    lineStrings <- mylines@Lines
    if(isTRUE(length(lineStrings) > 1)) {
      stop("Argument obj must have only one object of class Lines in each ",
           "position of the list in slot lines. Use ",
           "SpatialLinesToWKBMultiLineString instead of ",
           "SpatialLinesToWKBLineString.")
    }
    myline <- lineStrings[[1]]
    coords <- myline@coords
    writeBin(nrow(coords), rc, size = 4, endian = endian)
    apply(X = coords, MARGIN = 1, FUN = function(coord) {
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


#' Envelope of SpatialLines
#'
#' Takes an object of class \code{SpatialLines} or \code{SpatialLinesDataFrame}
#' and returns a data frame with six columns representing the envelope of each
#' object of class \code{Lines}.
#'
#' This function is called by the \code{\link{writeEnvelope}} function. Call the
#' \code{\link{writeEnvelope}} function instead of calling this function
#' directly.
#'
#' @param obj an object of class
#'   \code{\link[sp:SpatialLines-class]{SpatialLines}} or
#'   \code{\link[sp:SpatialLinesDataFrame-class]{SpatialLinesDataFrame}}.
#' @param centerfun function to apply to the x-axis limits and y-axis limits of
#'   the bounding box to obtain the x-coordinate and y-coordinate of the center
#'   of the bounding box.
#' @return A data frame with six columns named XMax, XMin, YMax, YMin, XCenter,
#'   and YCenter. The first four columns represent the corners of the bounding
#'   box of each object of class \code{Lines}. The last two columns represent
#'   the center of the bounding box of each object of class \code{Lines}. The
#'   number of rows in the returned data frame is the same as the length of the
#'   argument \code{obj}.
#'
#'   When this function is run in TIBCO Enterprise Runtime for R
#'   (\acronym{TERR}), the columns of the returned data frame have the
#'   SpotfireColumnMetaData attribute set to enable TIBCO Spotfire to recognize
#'   them as containing envelope information.
#' @seealso \code{\link{writeEnvelope}}
#'
#'   Example usage at \code{\link{SpatialLinesToWKBMultiLineString}}
#' @noRd
#' @importFrom sp bbox
SpatialLinesEnvelope <- function(obj, centerfun = mean) {
  if(is.character(centerfun)) {
    centerfun <- eval(parse(text = centerfun))
  }
  coords <- as.data.frame(t(vapply(X = obj@lines, FUN = function(mylines) {
    c(XMax = bbox(mylines)["x", "max"],
      XMin = bbox(mylines)["x", "min"],
      YMax = bbox(mylines)["y", "max"],
      YMin = bbox(mylines)["y", "min"],
      XCenter = centerfun(bbox(mylines)["x", ], na.rm = TRUE),
      YCenter = centerfun(bbox(mylines)["y", ], na.rm = TRUE))
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
