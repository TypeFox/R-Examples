#' Maximum of Absolute Values of Vector
#'
#' Take the maximum of the absolute values of an input vector.
#'
#' @param x A numeric vector or array.
#' @param na.rm A logical indicating whether missing values should be removed.
#'
#' @return The maximum of the absolute values of elements of the input vector.
#'
#' @export
#'
#' @examples
#' x <- c(5, 3, -9, -100, 3.14159, 7.5)
#' absmax(x)

absmax <- function(x, na.rm = FALSE) {
    max(abs(x), na.rm = na.rm)
}
