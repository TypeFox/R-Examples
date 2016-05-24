#' Make WKT multipoint objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @family R-objects
#' @examples
#' ## empty multipoint
#' multipoint("empty")
#' # multipoint("stuff")
#'
#' # numeric
#' multipoint(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180))
#'
#' # data.frame
#' df <- us_cities[1:25, c('long', 'lat')]
#' multipoint(df)
#'
#' # matrix
#' mat <- matrix(c(df$long, df$lat), ncol = 2)
#' multipoint(mat)
#'
#' # list
#' multipoint(list(c(100.000, 3.101), c(101.000, 2.100), c(3.140, 2.180)))
multipoint <- function(..., fmt = 16) {
  UseMethod("multipoint")
}

#' @export
multipoint.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('MULTIPOINT EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
multipoint.numeric <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  invisible(lapply(pts, checker, type = 'MULTIPOINT', len = 2))
  str <- paste0(lapply(pts, function(z){
    sprintf("(%s)", paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " "))
  }), collapse = ", ")
  sprintf('MULTIPOINT (%s)', str)
}

#' @export
multipoint.data.frame <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- paste0(apply(pts[[1]], 1, function(z){
    sprintf("(%s)", paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " "))
  }), collapse = ", ")
  sprintf('MULTIPOINT (%s)', str)
}

#' @export
multipoint.matrix <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- paste0(apply(pts[[1]], 1, function(z){
    sprintf("(%s)", paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " "))
  }), collapse = ", ")
  sprintf('MULTIPOINT (%s)', str)
}

#' @export
multipoint.list <- function(..., fmt = 16) {
  pts <- list(...)[[1]]
  fmtcheck(fmt)
  str <- paste0(lapply(pts, function(z) {
    sprintf("(%s)", paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " "))
  }), collapse = ", ")
  sprintf('MULTIPOINT (%s)', str)
}
