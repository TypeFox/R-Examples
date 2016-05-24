#' Make WKT linestring objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @family R-objects
#' @examples
#' ## empty linestring
#' linestring("empty")
#' # linestring("stuff")
#'
#' ## character string
#' linestring("LINESTRING (-116.4 45.2, -118.0 47.0)")
#'
#' # numeric
#' ## 2D
#' linestring(c(100.000, 0.000), c(101.000, 1.000), fmt=2)
#' linestring(c(100.0, 0.0), c(101.0, 1.0), c(120.0, 5.00), fmt=2)
#' ## 3D
#' linestring(c(0.0, 0.0, 10.0), c(2.0, 1.0, 20.0),
#'            c(4.0, 2.0, 30.0), c(5.0, 4.0, 40.0), fmt=2)
#' ## 4D
#' linestring(c(0.0, 0.0, 10.0, 5.0), c(2.0, 1.0, 20.0, 5.0),
#'            c(4.0, 2.0, 30.0, 5.0), c(5.0, 4.0, 40.0, 5.0), fmt=2)
#'
#' # data.frame
#' df <- data.frame(lon=c(-116.4,-118), lat=c(45.2,47))
#' linestring(df, fmt=1)
#' df <- data.frame(lon=c(-116.4,-118,-120), lat=c(45.2,47,49))
#' linestring(df, fmt=1)
#'
#' # matrix
#' mat <- matrix(c(-116.4,-118, 45.2, 47), ncol = 2)
#' linestring(mat, fmt=1)
#' mat2 <- matrix(c(-116.4, -118, -120, 45.2, 47, 49), ncol = 2)
#' linestring(mat2, fmt=1)
#'
#' # list
#' linestring(list(c(100.000, 0.000), c(101.000, 1.000)), fmt=2)
linestring <- function(..., fmt = 16) {
  UseMethod("linestring")
}

#' @export
linestring.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('LINESTRING EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
linestring.numeric <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  invisible(lapply(pts, checker, type = 'LINESTRING', len = 2:4))
  str <- paste0(lapply(pts, function(z){
    paste0(gsub("\\s", "", format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
  }), collapse = ", ")
  sprintf('LINESTRING (%s)', str)
}

#' @export
linestring.data.frame <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  str <- paste0(apply(pts[[1]], 1, function(x) paste0(format(x, nsmall = fmt, trim = TRUE), collapse = " ")), collapse = ", ")
  sprintf('LINESTRING (%s)', str)
}

#' @export
linestring.matrix <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  str <- paste0(apply(pts[[1]], 1, function(x) paste0(format(x, nsmall = fmt, trim = TRUE), collapse = " ")), collapse = ", ")
  sprintf('LINESTRING (%s)', str)
}

#' @export
linestring.list <- function(..., fmt = 16) {
  pts <- list(...)[[1]]
  fmtcheck(fmt)
  str <- paste0(lapply(pts, function(z){
    paste0(gsub("\\s", "", format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
  }), collapse = ", ")
  sprintf('LINESTRING (%s)', str)
}
