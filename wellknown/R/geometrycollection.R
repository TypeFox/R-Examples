#' Make WKT geometrycollection objects
#'
#' @export
#'
#' @param ... Character string WKT objects representing a Point, LineString,
#' Polygon, etc.
#' @details This is different from the other functions that create WKT from R
#' objects, in that we can't do the same thing for GeometryCollection's since
#' many different WkT object could be created from the same input. So,
#' this function accepts WKT strings already formed and attempts to creat a c
#' GeommetryCollection from them.
#' @family R-objects
#' @examples
#' ## empty geometrycollection
#' geometrycollection("empty")
#' # geometrycollection("stuff")
#'
#' # Character string, returns itself
#' geometrycollection("GEOMETRYCOLLECTION(POINT(4 6), LINESTRING(4 6, 7 10))")
#'
#' # From a point
#' geometrycollection(point(-116.4, 45.2))
#'
#' # From two points
#' geometrycollection(point(-116.4, 45.2), point(-118.4, 49.2))
#'
#' # From various object types
#' geometrycollection(point(-116.4, 45.2),
#'  linestring("LINESTRING (-116.4 45.2, -118.0 47.0)"),
#'  circularstring(list(c(1, 5), c(6, 2), c(7, 3)), fmt = 2)
#' )
geometrycollection <- function(...) {
  UseMethod("geometrycollection")
}

#' @export
geometrycollection.character <- function(...) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('GEOMETRYCOLLECTION EMPTY')
  } else {
    if (grepl("geometrycollection", pts[[1]], ignore.case = TRUE)) {
      check_str(pts[[1]])
    } else {
      if (!all(vapply(pts, lint, logical(1)))) stop("All inputs must be WKT strings", call. = FALSE)
      sprint("GEOMETRYCOLLECTION", paste(pts, collapse = ", "))
    }
  }
}
