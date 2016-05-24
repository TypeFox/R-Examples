#' @include source.R
#' @include osm-descriptors.R
#' @include as-osmar.R
{}



#' Get OSM data
#'
#' Get OSM data as \code{\link{osmar}} object from different sources
#' by providing a bounding box.
#'
#' @param x Data identifier, e.g., bounding box or specific element;
#'   see the help page of the used OSM source for a detailed list on
#'   the supported identifiers
#' @param source OSM source, e.g., \code{\link{osmsource_api}}
#' @param ... Additional arguments suppported by the specific OSM
#'   source; see corresponding source help page for a detailed list
#'
#' @return An \code{\link{osmar}} object
#'
#' @examples
#'   \dontrun{
#'   api <- osmsource_api()
#'
#'   box <- corner_bbox(11.579341, 48.15102, 11.582852, 48.1530)
#'   gschw <- get_osm(box, source = api)
#'
#'   kaufstr <- get_osm(way(3810479))
#'   kaufstr_full <- get_osm(way(3810479), full = TRUE)
#'   }
#'
#' @seealso \code{\link{bbox}}, \code{\link{osm_descriptors}},
#'   \code{\link{osmsource_api}}, \code{\link{osmsource_osmosis}}
#'
#' @import RCurl 
#' @import XML
#' @import gtools
#' @import methods
#' 
#' @export
get_osm <- function(x, source = osmsource_api(), ...) {
  raw <- get_osm_data(source, x, ...)

  ret <- xmlParse(raw)
  ret <- as_osmar(ret)
  #attr(ret, "identifier") <- x
  #attr(ret, "source") <- source

  ret
}



### Bounding box: ####################################################


#' Get OSM elements
#'
#' Utility functions to specify \emph{what} to get from the OSM data
#' source. These are the request elements which work for most sources,
#' see the specific sources for specialized elements.
#'
#' @param left Minimum longitude
#' @param bottom Minimum latitude
#' @param right Maximum longitude
#' @param top Maximum latitutde
#'
#' @seealso \code{\link{osm_descriptors}}, \code{\link{get_osm}}
#'
#' @aliases bbox
#' @rdname bbox
#' @family as_osmar_bbox
#'
#' @export
corner_bbox <- function(left, bottom, right, top) {
  ## TODO: check arguments
  structure(c(left = left, bottom = bottom,
              right = right, top = top), class = "bbox")
}



#' @param center_lon Center longitude
#' @param center_lat Center latitude
#' @param width Box width
#' @param height Box height
#'
#' @rdname bbox
#'
#' @export
center_bbox <- function(center_lon, center_lat, width, height) {
  stopifnot(center_lon <= 180 & center_lon >= -180)
  stopifnot(center_lat <= 90 & center_lat >= -90)

  width <- width / 2
  height <- height / 2

  a <- 6378137
  esq <- (2 - (1/298.257223563)) * (1/298.257223563)
  W <- sqrt(1 - esq * (sin(center_lat * pi/180))^2)
  M <- a * (1 - esq)/W^3
  mPerLatD <- 1/((pi/180) * M)
  top <- center_lat + mPerLatD * height
  bottom <- center_lat - mPerLatD * height
  N <- a/W
  mPerLonD <- 1/((pi/180) * N * cos(center_lat * pi/180))
  left <- center_lon - mPerLonD * width
  right <- center_lon + mPerLonD * width

  if (left < -180)
    left <- left + 360
  if (right > 180)
    right <- right - 360

  corner_bbox(left, bottom, right, top)
}



#' Bounding box converter generic
#'
#' Generic function for implementing converters from various objects
#' (e.g., \link[sp]{sp} \code{\link[sp]{Spatial}} objects) to osmar
#' \code{\link{bbox}} objects.
#'
#' @param obj Object to compute osmar \code{\link{bbox}}
#' @param ... Additional parameters for underlying functions
#'
#' @family as_osmar_bbox
#'
#' @export
as_osmar_bbox <- function(obj, ...) {
  UseMethod("as_osmar_bbox")
}



size <- function(x, ...) {
  UseMethod("size")
}

size.bbox <- function(x) {
  unname((x[1] - x[3]) * (x[2] - x[4]))
}


