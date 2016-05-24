#' Make WKT multipolygon objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @details There is no \code{numeric} input option for multipolygon. There is no way as of
#' yet to make a nested multipolygon with \code{data.frame} input, but you can do so
#' with list input. See examples.
#' @family R-objects
#' @examples
#' ## empty multipolygon
#' multipolygon("empty")
#' # multipolygon("stuff")
#'
#' # data.frame
#' df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20))
#' df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5))
#' multipolygon(df, df2, fmt=0)
#' multipolygon(df, df2, fmt=0) %>% lint
#' multipolygon(df, df2) %>% wktview(zoom=3)
#'
#' # matrix
#' mat <- matrix(c(df$long, df$lat), ncol = 2)
#' mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
#' multipolygon(mat, mat2, fmt=0)
#'
#' # list
#' multipolygon(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20)),
#'   list(c(15, 5), c(40, 10), c(10, 20), c(5, 10), c(15, 5)), fmt=2)
#'
#' polys <- list(
#'   list(c(30, 20), c(45, 40), c(10, 40), c(30, 20)),
#'   list(c(15, 5), c(40, 10), c(10, 20), c(5, 10), c(15, 5))
#' )
#' multipolygon(polys, fmt=2) %>%
#'   wktview(zoom=3)
#'
#' ## nested polygons
#' polys <- list(
#'   list(c(40, 40), c(20, 45), c(45, 30), c(40, 40)),
#'   list(
#'     list(c(20, 35), c(10, 30), c(10, 10), c(30, 5), c(45, 20), c(20, 35)),
#'     list(c(30, 20), c(20, 15), c(20, 25), c(30, 20))
#'   )
#' )
#' multipolygon(polys, fmt=0)
#' multipolygon(polys, fmt=0) %>% lint
multipolygon <- function(..., fmt = 16) {
  UseMethod("multipolygon")
}

#' @export
multipolygon.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('MULTIPOLYGON EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
multipolygon.data.frame <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- lapply(pts, function(v) {
    sprintf("((%s))", paste0(apply(v, 1, function(z){
      paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
    }), collapse = ", "))
  })
#   str <- paste0(apply(pts[[1]], 1, function(z){
#     paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
#   }), collapse = ", ")
  sprintf('MULTIPOLYGON (%s)', paste0(str, collapse = ", "))
}

#' @export
multipolygon.matrix <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- lapply(pts, function(v) {
    sprintf("((%s))", paste0(apply(v, 1, function(z){
      paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
    }), collapse = ", "))
  })
  #   str <- paste0(apply(pts[[1]], 1, function(z){
  #     paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
  #   }), collapse = ", ")
  sprintf('MULTIPOLYGON (%s)', paste0(str, collapse = ", "))
}

#' @export
multipolygon.list <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  pts <- un_nest(pts)
  str <- lapply(pts, function(z) {
    if (length(z) > 1 && sapply(z, class)[1] != "numeric") {
      inparens(paste0(lapply(z, make1multipoly, fmt = fmt), collapse = ", "))
    } else {
      inparens(make1multipoly(z, fmt))
    }
  })
  sprintf('MULTIPOLYGON (%s)', paste0(str, collapse = ", "))
}

make1multipoly <- function(m, fmt) {
  inparens(paste0(lapply(m, function(b) paste0(str_trim_(format(b, nsmall = fmt, trim = TRUE)), collapse = " ")), collapse = ", "))
}

inparens <- function(x) {
  sprintf("(%s)", x)
}
