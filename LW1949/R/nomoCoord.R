#' Find the Coordinate from the Scale of a Nomograph
#'
#' Find the x and y coordinate corresponding to a specified point on a
#' nomograph scale.
#' @param df
#'   A data frame with three columns, the \code{x} and \code{y} coordinates
#'   of the nomograph scale and the corresponding values.
#' @param val
#'   A numeric scalar identifying the point on the scale for which the
#'   coordinates will be returned.
#' @export
#' @import
#'   plotrix
#' @details
#' The function makes it easier to add points or lines to a nomograph for
#' illustrative purposes.
#'
#' Each scale is assumed to be displayed on the log10 scale.
#' @return
#'   A numeric vector of length two with the x and y coordinates of the
#'   specified point in plotting units of the nomograph.
#' @examples
#' scales <- LWnomo1(TRUE)
#' fromxy <- nomoCoord(scales$scale1r, 34)
#' toxy <- nomoCoord(scales$scale3, 16^2/(100*34))
#' segments(fromxy[1], fromxy[2], toxy[1], toxy[2], col="red")
#'
nomoCoord <- function(df, val) {
  if(val < df$values[1] | val > df$values[2]) {
    stop("val must be between ", df$values[1], " and  ", df$values[1], ".")
  }
  c(df$x[1], rescale(log10(c(val, df$values)), df$y)[1])
}
