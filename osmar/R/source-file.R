#' @include get.R
#' @include source.R
{}



#' OSM file data source
#'
#' Imports the complete OSM file.
#'
#' @section Supported request elements:
#'
#' \describe{
#'
#'   \item{Dummy request element:}{Use the function
#'     \code{compete_file} as dummy description for all elements
#      available in the file.
#'   }
#'
#' }
#'
#' @param file The file name (and path) of the osm file
#'
#' @examples
#'   \dontrun{
#'     get_osm(complete_file(), source = osmsource_file("muc.osm"))
#'   }
#'
#' @seealso \code{\link{get_osm}}, \code{\link{bbox}},
#'   \code{\link{osm_descriptors}}
#' @aliases complete_file
#' @family osmsource
#'
#' @export
osmsource_file <- function(file) {
  osmsource(list(file = file), "osmfile")
}



get_osm_data.osmfile <- function(source, what, ...) {
  readLines(source$file)
}



#' @export
complete_file <- function() {
  structure(c(from = 0, to = Inf), class = "complete_file")
}
