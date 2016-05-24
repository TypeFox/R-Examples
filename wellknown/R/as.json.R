#' Convert geojson R list to JSON
#'
#' @export
#' @importFrom jsonlite toJSON
#'
#' @param x Output from \code{\link{wkt2geojson}}
#' @param pretty (logical) Adds indentation whitespace to JSON output. Can be TRUE/FALSE
#' or a number specifying the number of spaces to indent. See \code{\link[jsonlite]{prettify}}
#' Default: \code{TRUE} Having \code{TRUE} as default makes it easy to copy paste to a
#' text editor, etc.
#' @param auto_unbox (logical) Automatically unbox all atomic vectors of length 1.
#' Default: \code{TRUE}
#' @param ... Further args passed on to \code{\link[jsonlite]{toJSON}}
#' @examples
#' str <- "POLYGON ((100 0.1, 101.1 0.3, 101 0.5, 100 0.1),
#'    (103.2 0.2, 104.8 0.2, 100.8 0.8, 103.2 0.2))"
#' as_json(wkt2geojson(str))
#' as_json(wkt2geojson(str), FALSE)
as_json <- function(x, pretty=TRUE, auto_unbox=TRUE, ...) {
  UseMethod("as_json")
}

#' @export
as_json.geojson <- function(x, pretty = TRUE, auto_unbox = TRUE, ...) {
  jsonlite::toJSON(unclass(x), pretty = pretty, auto_unbox = auto_unbox, ...)
}
