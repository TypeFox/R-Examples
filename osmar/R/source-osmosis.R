#' @include get.R
#' @include source.R
{}



#' Osmosis OSM data source
#'
#' Planet dumps as OSM data source through the osmosis command line
#' Java application.
#'
#' Osmosis is a command line Java application for processing OSM
#' data. It allows, among other things, to extract data inside a
#' bounding box or polygon from so called planet dumps. The usage of
#' this source requires an installed osmosis; see
#' \url{http://wiki.openstreetmap.org/wiki/Osmosis}.
#'
#' @section Supported request elements:
#'
#' \describe{
#'
#'   \item{Bounding box:}{
#'
#'     Use \code{\link{corner_bbox}} or \code{\link{center_bbox}} to
#'     retrieve:
#'
#'     \itemize{
#'
#'       \item all nodes that are inside a given bounding box and any
#'         relations that reference them;
#'
#'       \item all ways that reference at least one node that is
#'         inside a given bounding box, any relations that reference
#'         them [the ways], and any nodes outside the bounding box
#'         that the ways may reference;
#'
#'       \item all relations that reference one of the nodes or ways
#'         included due to the above rules (does not apply
#'         recursively);
#'     }
#'
#'  }
#'
#' }
#'
#' @param file The file name (and path) of the planet dump
#' @param osmosis The path to the osmosis application
#'
#' @examples
#'   \dontrun{
#'     ## Download and extract a planet file:
#'     download.file("http://osmar.r-forge.r-project.org/",
#'                   "muenchen.osm.gz")
#'     system("gzip -d muenchen.osm.gz")
#'
#'     ## Define osmosis source; note that we assume that
#'     ## osmosis is in our path environment variable (if
#'     ## not, set osmosis argument to the executable):
#'     src <- osmsource_osmosis(file = "muenchen.osm")
#'
#'     ## Get the center of Munich:
#'     muc_bbox <- center_bbox(11.575278, 48.137222,
#'                             3000, 3000)
#'     muc <- get_osm(muc_bbox, src)
#'     muc
#'   }
#'
#' @references
#'   \url{http://wiki.openstreetmap.org/wiki/Osmosis}
#' @seealso \code{\link{get_osm}}, \code{\link{bbox}},
#'   \code{\link{osm_descriptors}}
#' @family osmsource
#'
#' @export
osmsource_osmosis <- function(file, osmosis = "osmosis") {
  osmsource(list(file = file, osmosis = osmosis), "osmosis")
}



get_osm_data.osmosis <- function(source, what, ...) {
  destination <- tempfile()

  request <- osm_request(source, what, destination)
  ret <- system(request, ...)
  response <- readLines(destination)

  unlink(destination)

  response
}



setMethod("osm_request", signature = c("osmosis", "bbox"),
function(source, what, destination, ...) {
  fin <- sprintf("--read-xml enableDateParsing=no file=%s", source$file)
  fout <- sprintf("--write-xml file=%s", destination)

  args <- sprintf("--bounding-box top=%s left=%s bottom=%s right=%s",
                  what["top"], what["left"], what["bottom"], what["right"])

  sprintf("%s %s %s %s", source$osmosis, fin, args, fout)
})



setMethod("osm_request", signature = c("osmosis", "element"),
function(source, what, destination, ...) {
  stop("Not implemented yet.")

  #fin <- sprintf("--read-xml enableDateParsing=no file=%s", source$file)
  #fout <- sprintf("--write-xml file=%s", destination)

  #args <- ""

  #sprintf("%s %s %s %s", source$osmosis, fin, args, fout)
})




### Special omsosis requester: #######################################

setOldClass(c("osmosis_args"))

osmosis_args <- function(...) {
  structure(paste(..., collapse = " "), class = "osmosis_args")
}

setMethod("osm_request", signature = c("osmosis", "osmosis_args"),
function(source, what, destination, ...) {
  fin <- sprintf("--read-xml enableDateParsing=no file=%s", source$file)
  fout <- sprintf("--write-xml file=%s", destination)

  args <- unclass(what)

  sprintf("%s %s %s %s", source$osmosis, fin, args, fout)
})


