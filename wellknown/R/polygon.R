#' Make WKT polygon objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon, MultiPolygon, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @details You can create nested polygons with \code{list} and \code{data.frame} inputs,
#' but not from \code{numeric} inputs. See examples.
#' @family R-objects
#' @examples
#' ## empty polygon
#' polygon("empty")
#' # polygon("stuff")
#'
#' # numeric
#' polygon(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001), fmt=2)
#'
#' # data.frame
#' ## single polygon
#' df <- us_cities[2:5,c('long','lat')]
#' df <- rbind(df, df[1,])
#' polygon(df, fmt=2) %>% wktview(zoom=4)
#' ## nested polygons
#' df2 <- data.frame(long = c(-85.9, -85.9, -93, -93, -85.9),
#'                   lat = c(37.5, 35.3, 35.3, 37.5, 37.5))
#' polygon(df, df2, fmt=2) %>% wktview(zoom=4)
#'
#' # matrix
#' mat <- matrix(c(df$long, df$lat), ncol = 2)
#' polygon(mat)
#'
#' # list
#' # single list - creates single polygon
#' ply <- list(c(100.001, 0.001), c(101.12345, 0.001), c(101.001, 1.001), c(100.001, 0.001))
#' polygon(ply, fmt=2) %>% wktview(zoom=7)
#' # nested list - creates nested polygon
#' polygon(list(c(35, 10), c(45, 45), c(15, 40), c(10, 20), c(35, 10)),
#'    list(c(20, 30), c(35, 35), c(30, 20), c(20, 30)), fmt=2) %>%
#'    wktview(zoom=3)
#' # multiple lists nested within a list
#' polygon(list(list(c(35, 10), c(45, 45), c(15, 40), c(10, 20), c(35, 10)),
#'    list(c(20, 30), c(35, 35), c(30, 20), c(20, 30))), fmt=2) %>%
#'    wktview(zoom=3)
polygon <- function(..., fmt = 16) {
  UseMethod("polygon")
}

#' @export
polygon.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('POLYGON EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
polygon.numeric <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  invisible(lapply(pts, checker, type = 'POLYGON', len = 2))
  str <- paste0(lapply(pts, function(z){
    paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
  }), collapse = ", ")
  sprintf('POLYGON ((%s))', str)
}

#' @export
polygon.data.frame <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- lapply(pts, function(v) {
    sprintf("(%s)", paste0(apply(v, 1, function(z){
      paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
    }), collapse = ", "))
  })
  sprintf('POLYGON (%s)', paste0(str, collapse = ", "))
}

#' @export
polygon.matrix <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  # invisible(lapply(pts, checker, type='MULTIPOINT', len=2))
  str <- lapply(pts, function(v) {
    sprintf("(%s)", paste0(apply(v, 1, function(z){
      paste0(str_trim_(format(z, nsmall = fmt, trim = TRUE)), collapse = " ")
    }), collapse = ", "))
  })
  sprintf('POLYGON (%s)', paste0(str, collapse = ", "))
}

#' @export
polygon.list <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  pts <- un_nest(pts)
  str <- sprintf("(%s)", lapply(pts, function(z) {
    paste0(lapply(z, function(b) paste0(str_trim_(format(b, nsmall = fmt, trim = TRUE)), collapse = " ")), collapse = ", ")
  }))
  sprintf('POLYGON (%s)', paste0(str, collapse = ", "))
}

un_nest <- function(x) {
  first <- sapply(x, class)
  if (length(first) == 1 && first == "list") {
    if (sapply(x[[1]], class)[1] == "list") {
      unlist(x, recursive = FALSE)
    } else {
      return(x)
    }
  } else {
    return(x)
  }
}
