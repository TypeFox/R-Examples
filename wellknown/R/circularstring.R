#' Make WKT circularstring objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @family R-objects
#' @examples
#' ## empty circularstring
#' circularstring("empty")
#' # circularstring("stuff")
#'
#' # Character string
#' circularstring("CIRCULARSTRING(1 5, 6 2, 7 3)")
#'
#' # data.frame
#' df <- data.frame(lon = c(-116.4, -118), lat = c(45.2, 47))
#' circularstring(df, fmt=1)
#' df <- data.frame(lon=c(-116.4, -118, -120), lat=c(45.2, 47, 49))
#' circularstring(df, fmt=1)
#'
#' # matrix
#' mat <- matrix(c(-116.4,-118, 45.2, 47), ncol = 2)
#' circularstring(mat, fmt=1)
#' mat2 <- matrix(c(-116.4, -118, -120, 45.2, 47, 49), ncol = 2)
#' circularstring(mat2, fmt=1)
#'
#' # list
#' x <- list(c(1, 5), c(6, 2), c(7, 3))
#' circularstring(x, fmt=2)
circularstring <- function(..., fmt = 16) {
  UseMethod("circularstring")
}

#' @export
circularstring.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('CIRCULARSTRING EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
circularstring.data.frame <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  str <- p0c(apply(pts[[1]], 1, function(x) p0c(format(x, nsmall = fmt, trim = TRUE))), cl = ", ")
  sprint('CIRCULARSTRING', str)
}

#' @export
circularstring.matrix <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  str <- p0c(apply(pts[[1]], 1, function(x) p0c(format(x, nsmall = fmt, trim = TRUE))), cl = ", ")
  sprint('CIRCULARSTRING', str)
}

#' @export
circularstring.list <- function(..., fmt = 16) {
  pts <- list(...)[[1]]
  fmtcheck(fmt)
  str <- paste0(lapply(pts, function(z){
    paste0(gsub("\\s", "", format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
  }), collapse = ", ")
  sprint('CIRCULARSTRING', str)
}

sprint <- function(type, str) {
  sprintf('%s (%s)', type, str)
}

p0c <- function(..., cl = " ") {
  paste0(..., collapse = cl)
}
