# Convert a well-known binary (WKB) geometry representation to an R spatial
# object

#' Convert \acronym{WKB} to Spatial Objects
#'
#' Converts well-known binary (\acronym{WKB}) geometry representations to
#' \code{Spatial} objects.
#'
#' @param wkb \code{list} in which each element is a \code{\link[base]{raw}}
#'   vector consisting of a \acronym{WKB} geometry representation.
#' @param id character vector of unique identifiers of geometries. The length of
#'   \code{id} must be the same as the length of the \code{wkb} list.
#' @param proj4string projection string of class
#'   \code{\link[sp:CRS-class]{CRS}}.
#' @details Supported \acronym{WKB} geometry types are Point, LineString,
#'   Polygon, MultiPoint, MultiLineString, and MultiPolygon. All elements in the
#'   \code{list} must have the same \acronym{WKB} geometry type. The
#'   \acronym{WKB} geometry representations may use little-endian or big-endian
#'   byte order.
#'
#'   The argument \code{wkb} may also be a \code{\link[base]{raw}} vector
#'   consisting of one \acronym{WKB} geometry representation. In that case, the
#'   argument \code{id} must have length one.
#' @return An object inheriting class \code{\link[sp:Spatial-class]{Spatial}}.
#'
#'   The return value may be an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}},
#'   \code{\link[sp:SpatialLines-class]{SpatialLines}},
#'   \code{\link[sp:SpatialPolygons-class]{SpatialPolygons}}, or a \code{list}
#'   in which each element is an object of class
#'   \code{\link[sp:SpatialPoints-class]{SpatialPoints}}. The class of the
#'   return value depends on the \acronym{WKB} geometry type as shown in the
#'   table below.
#'
#'   \tabular{ll}{
#'   \strong{Type of \acronym{WKB} geometry} \tab \strong{Class of return value}\cr
#'   Point \tab \code{SpatialPoints}\cr
#'   LineString \tab \code{SpatialLines}\cr
#'   Polygon \tab \code{SpatialPolygons}\cr
#'   MultiPoint \tab \code{list} of \code{SpatialPoints}\cr
#'   MultiLineString \tab \code{SpatialLines}\cr
#'   MultiPolygon \tab \code{SpatialPolygons}\cr
#'   }
#' @examples
#' # create a list of WKB geometry representations of type Point
#' wkb <- list(
#'   as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
#'            0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40)),
#'   as.raw(c(0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
#'            0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40))
#' )
#'
#' # convert to object of class SpatialPoints
#' obj <- readWKB(wkb)
#'
#'
#' # create a list of WKB geometry representations of type MultiPoint
#' wkb <- list(
#'   as.raw(c(0x01, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x01,
#'            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f,
#'            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x40)),
#'   as.raw(c(0x01, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x01,
#'            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40,
#'            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40)))
#'
#' # convert to list of objects of class SpatialPoints
#' obj <- readWKB(wkb)
#'
#'
#' # create a list of WKB geometry representations of type MultiLineString
#' wkb <- list(
#'   as.raw(c(0x01, 0x05, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x02,
#'            0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
#'            0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x08,
#'            0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00,
#'            0x00, 0x00, 0x00, 0x00, 0x00, 0x40)),
#'   as.raw(c(0x01, 0x05, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x02,
#'            0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
#'            0x00, 0x00, 0xf0, 0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0,
#'            0x3f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00,
#'            0x00, 0x00, 0x00, 0x00, 0xf8, 0x3f)))
#'
#' # convert to object of class SpatialLines
#' obj <- readWKB(wkb)
#'
#'
#' # create a list of WKB geometry representations of type Polygon
#' wkb <- list(
#'   as.raw(c(0x01, 0x03, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x05, 0x00,
#'            0x00, 0x00, 0x34, 0x03, 0xf0, 0xac, 0xce, 0x66, 0x5d, 0xc0, 0x8f,
#'            0x27, 0x95, 0x21, 0xab, 0xa6, 0x44, 0x40, 0xa0, 0x32, 0x81, 0x18,
#'            0x78, 0x83, 0x5d, 0xc0, 0xc8, 0xd2, 0xa0, 0xee, 0x23, 0x0b, 0x41,
#'            0x40, 0x80, 0xec, 0x72, 0x54, 0xde, 0xb1, 0x5f, 0xc0, 0xc8, 0xd2,
#'            0xa0, 0xee, 0x23, 0x0b, 0x41, 0x40, 0xec, 0x1b, 0x04, 0xc0, 0x87,
#'            0xce, 0x5f, 0xc0, 0x8f, 0x27, 0x95, 0x21, 0xab, 0xa6, 0x44, 0x40,
#'            0x34, 0x03, 0xf0, 0xac, 0xce, 0x66, 0x5d, 0xc0, 0x8f, 0x27, 0x95,
#'            0x21, 0xab, 0xa6, 0x44, 0x40)),
#'   as.raw(c(0x01, 0x03, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x05, 0x00,
#'            0x00, 0x00, 0x08, 0x36, 0xdc, 0x8b, 0x9f, 0x3d, 0x51, 0xc0, 0x0f,
#'            0xb3, 0x2a, 0x6a, 0x3f, 0x1c, 0x46, 0x40, 0x47, 0xcb, 0x54, 0xe7,
#'            0xcb, 0x5e, 0x51, 0xc0, 0x45, 0x81, 0x50, 0x31, 0xfa, 0x80, 0x42,
#'            0x40, 0xa9, 0xba, 0x74, 0x6d, 0xf5, 0xa1, 0x53, 0xc0, 0x45, 0x81,
#'            0x50, 0x31, 0xfa, 0x80, 0x42, 0x40, 0xe8, 0x4f, 0xed, 0xc8, 0x21,
#'            0xc3, 0x53, 0xc0, 0x0f, 0xb3, 0x2a, 0x6a, 0x3f, 0x1c, 0x46, 0x40,
#'            0x08, 0x36, 0xdc, 0x8b, 0x9f, 0x3d, 0x51, 0xc0, 0x0f, 0xb3, 0x2a,
#'            0x6a, 0x3f, 0x1c, 0x46, 0x40)))
#'
#' # convert to object of class SpatialPolygons
#' obj <- readWKB(wkb)
#'
#'
#' # specify id and proj4string
#' obj <- readWKB(
#'   wkb,
#'   id = c("San Francisco", "New York"),
#'   proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#' )
#' @seealso \code{\link{writeWKB}}, \code{\link{hex2raw}}
#' @keywords wkb
#' @export
#' @importFrom sp CRS
#' @importFrom sp SpatialPoints
#' @importFrom sp SpatialLines
#' @importFrom sp SpatialPolygons
readWKB <- function(wkb, id = NULL, proj4string = CRS(as.character(NA))) {
  if(inherits(wkb, "raw") && (is.null(id) || length(id) == 1)) {
    wkb <- list(wkb)
  }
  if(!is.list(wkb)) {
    stop("wkb must be a list")
  }
  if(isTRUE(length(wkb) < 1)) {
    stop("wkb must have length 1 or greater")
  }
  if(!all(vapply(X = wkb, FUN = inherits, FUN.VALUE = logical(1), "raw"))) {
    stop("Each element of wkb must be a raw vector")
  }
  if(is.null(id)) {
    id <- as.character(seq_along(wkb))
  }
  if(!identical(length(wkb), length(id))) {
    stop("wkb and id must have same length")
  }
  if(is.character(proj4string)) {
    proj4string = CRS(proj4string)
  }
  if(length(proj4string) != 1) {
    stop("proj4string must have length 1")
  }
  obj <- mapply(wkb, id, FUN = function(WkbGeom, Id) {
    rc <- rawConnection(WkbGeom, "r")
    on.exit(close(rc))
    seek(rc, 0L)
    byteOrder <- readByteOrder(rc)
    if(byteOrder == as.raw(1L)) {
      endian <- "little"
    } else {
      endian <- "big"
    }
    wkbType <- readWkbType(rc, endian)
    if(wkbType == 1L) {
      readWkbPoint(rc, Id, endian)
    } else if(wkbType == 2L) {
      readWkbLineString(rc, Id, endian)
    } else if(wkbType == 3L) {
      readWkbPolygon(rc, Id, endian)
    } else if(wkbType == 4L) {
      readWkbMultiPoint(rc, Id, endian)
    } else if(wkbType == 5L) {
      readWkbMultiLineString(rc, Id, endian)
    } else if(wkbType == 6L) {
      readWkbMultiPolygon(rc, Id, endian)
    } else if(wkbType == 7L) {
      stop("GeometryCollection is not a supported geometry type")
    } else {
      stop("Supported geometry types are Point, LineString, Polygon, MultiPoint, MultiLineString, and MultiPolygon")
    }
  }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  objClass <- unique(vapply(obj, class, character(1)))
  if(isTRUE(length(objClass) > 1)) {
    stop("Elements of wkb cannot have different geometry types")
  }
  if(objClass == "numeric") {
    SpatialPoints(do.call("rbind", obj), proj4string = proj4string)
  } else if(objClass == "matrix" || objClass == "data.frame") {
    lapply(X = obj, FUN = SpatialPoints, proj4string = proj4string)
  } else if(objClass == "Lines") {
    SpatialLines(obj, proj4string = proj4string)
  } else if(objClass == "Polygons") {
    SpatialPolygons(obj, proj4string = proj4string)
  } else {
    stop("Unexpected object")
  }
}

readWkbMultiPoint <- function(rc, multiPointId, endian) {
  numPoints <- readInteger(rc, endian)
  t(vapply(seq_len(numPoints), function(...) {
    byteOrder <- readByteOrder(rc)
    if(byteOrder == as.raw(1L)) {
      endian <- "little"
    } else {
      endian <- "big"
    }
    wkbType <- readWkbType(rc, endian)
    if(wkbType != 1L) {
      stop("MultiPoints may contain only Points")
    }
    readPoint(rc, endian)
  }, numeric(2)))
}

#' @importFrom sp Lines
#' @importFrom sp Line
readWkbMultiLineString <- function(rc, multiLineStringId, endian) {
  numLineStrings <- readInteger(rc, endian)
  Lines(unlist(lapply(seq_len(numLineStrings), function(...) {
    byteOrder <- readByteOrder(rc)
    if(byteOrder == as.raw(1L)) {
      endian <- "little"
    } else {
      endian <- "big"
    }
    wkbType <- readWkbType(rc, endian)
    if(wkbType != 2L) {
      stop("MultiLineStrings may contain only LineStrings")
    }
    numPoints <- readInteger(rc, endian)
    Line(readPoints(rc, numPoints, endian))
  })), multiLineStringId)
}

#' @importFrom sp Polygons
readWkbMultiPolygon <- function(rc, multiPolygonId, endian) {
  numPolygons <- readInteger(rc, endian)
  Polygons(unlist(lapply(seq_len(numPolygons), function(...) {
    byteOrder <- readByteOrder(rc)
    if(byteOrder == as.raw(1L)) {
      endian <- "little"
    } else {
      endian <- "big"
    }
    wkbType <- readWkbType(rc, endian)
    if(wkbType != 3L) {
      stop("MultiPolygons may contain only Polygons")
    }
    numRings <- readInteger(rc, endian)
    readLinearRings(rc, numRings, endian)
  })), multiPolygonId)
}

readWkbPoint <- function(rc, pointId, endian) {
  readPoint(rc, endian)
}

#' @importFrom sp Lines
#' @importFrom sp Line
readWkbLineString <- function(rc, lineStringId, endian) {
  numPoints <- readInteger(rc, endian)
  Lines(list(Line(readPoints(rc, numPoints, endian))), lineStringId)
}

#' @importFrom sp Polygons
readWkbPolygon <- function(rc, polygonId, endian) {
  numRings <- readInteger(rc, endian)
  Polygons(readLinearRings(rc, numRings, endian), polygonId)
}

readLinearRings <- function(rc, numRings, endian) {
  lapply(seq_len(numRings), function(ringId) {
    readLinearRing(rc, endian)
  })
}

#' @importFrom sp Polygon
readLinearRing <- function(rc, endian) {
  numPoints <- readInteger(rc, endian)
  Polygon(readPoints(rc, numPoints, endian))
}

readPoints <- function(rc, numPoints, endian) {
  t(vapply(seq_len(numPoints), function(pointId) {
    readPoint(rc, endian)
  }, numeric(2)))
}

readPoint <- function(rc, endian) {
  c(x = readDouble(rc, endian), y = readDouble(rc, endian))
}

readByteOrder <- function(rc) {
  readByte(rc)
}

readWkbType <- function(rc, endian) {
  readInteger(rc, endian)
}

readByte <- function(rc) {
  readBin(rc, what = "raw", size = 1L)
}

readInteger <- function(rc, endian) {
  readBin(rc, what = "integer", size = 4L, endian = endian)
}

readDouble <- function(rc, endian) {
  readBin(rc, what = "double", size = 8L, endian = endian)
}
