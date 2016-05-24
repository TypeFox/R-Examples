#' Make WKT multilinestring objects
#'
#' @export
#'
#' @param ... A GeoJSON-like object representing a Point, LineString, Polygon,
#' multilinestring, etc.
#' @param fmt Format string which indicates the number of digits to display after the
#' decimal point when formatting coordinates. Max: 20
#' @details There is no \code{numeric} input option for multilinestring. There is no
#' way as of yet to make a nested multilinestring with \code{data.frame} input, but you
#' can do so with list input. See examples.
#' @family R-objects
#' @examples
#' ## empty multilinestring
#' multilinestring("empty")
#' # multilinestring("stuff")
#'
#' # character string
#' x <- "MULTILINESTRING ((30 20, 45 40, 10 40), (15 5, 40 10, 10 20))"
#' multilinestring(x)
#'
#' # data.frame
#' df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40))
#' df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20))
#' multilinestring(df, df2, fmt=0)
#' multilinestring(df, df2, fmt=0) %>% lint
#' multilinestring(df, df2) %>% wktview(zoom=3)
#'
#' # matrix
#' mat <- matrix(c(df$long, df$lat), ncol = 2)
#' mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
#' multilinestring(mat, mat2, fmt=0)
#'
#' # list
#' x1 <- list(c(30, 20), c(45, 40), c(10, 40))
#' x2 <- list(c(15, 5), c(40, 10), c(10, 20))
#' multilinestring(x1, x2, fmt=2)
#'
#' polys <- list(
#'   list(c(30, 20), c(45, 40), c(10, 40)),
#'   list(c(15, 5), c(40, 10), c(10, 20))
#' )
#' multilinestring(polys, fmt=2) %>%
#'   wktview(zoom=3)
multilinestring <- function(..., fmt = 16) {
  UseMethod("multilinestring")
}

#' @export
multilinestring.character <- function(..., fmt = 16) {
  pts <- list(...)
  if (grepl("empty", pts[[1]], ignore.case = TRUE)) {
    return('MULTILINESTRING EMPTY')
  } else {
    check_str(pts)
  }
}

#' @export
multilinestring.data.frame <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  str <- lapply(pts, function(v) {
    sprintf("(%s)", paste0(apply(v, 1, function(z){
      p0c(str_trim_(format(z, nsmall = fmt, trim = TRUE)))
    }), collapse = ", "))
  })
  sprint_multil(str)
}

#' @export
multilinestring.matrix <- function(..., fmt = 16){
  pts <- list(...)
  fmtcheck(fmt)
  str <- lapply(pts, function(v) {
    sprintf("(%s)", paste0(apply(v, 1, function(z){
      p0c(str_trim_(format(z, nsmall = fmt, trim = TRUE)))
    }), collapse = ", "))
  })
  sprint_multil(str)
}

#' @export
multilinestring.list <- function(..., fmt = 16) {
  pts <- list(...)
  fmtcheck(fmt)
  pts <- un_nest(pts)
  str <- lapply(pts, function(z) {
    if (length(z) > 1 && sapply(z, class)[1] != "numeric") {
      inparens(paste0(lapply(z, make1multipoly, fmt = fmt), collapse = ", "))
    } else {
      make1multilinestr(z, fmt)
    }
  })
  sprint_multil(str)
}

sprint_multil <- function(x) {
  sprintf('MULTILINESTRING (%s)', paste0(x, collapse = ", "))
}

sprint <- function(type, str) {
  sprintf('%s (%s)', type, str)
}

make1multilinestr <- function(m, fmt) {
  inparens(paste0(lapply(m, function(b) {
      paste0(str_trim_(format(b, nsmall = fmt, trim = TRUE)), collapse = " ")
    }), collapse = ", ")
  )
}
