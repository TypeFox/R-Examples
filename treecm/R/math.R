#' @title Computes the x,y cartesian coordinates
#'
#' @description Computes the \eqn{x} and \eqn{y} cartesian coordinates from a set of polar coordinates
#'
#' @param angle The angle in degrees (measured clockwise from the North or  
#' any other relevant bearing system defined in the field)
#' @param distance The distance
#' @return A vector holding the \eqn{x} and \eqn{y} coordinats expressed in the same unit as the distance argument
#'
#' @note The function assumes the angle is measured clockwise whereas
#' trigonometric functions require a conventional counterclockwise 
#' measured angle. Thus the function computes \eqn{x} coordinate as the sine of 
#' the angle, and the \eqn{y} coordinate as the cosine of 
#' the angle, enabling a correct representation of them on a cartesian plot.
#' @export
toCartesianXY <- function(angle, distance) {
  angleRad <- angle * pi / 180
  c(sin(angleRad) * distance, cos(angleRad) * distance)
}


#' @title Converts cartesian (x, y) into polar (angle, distance) coordinates
#'
#' @description Converts cartesian coordinates (\eqn{x}, \eqn{y}) into polar (angle, distance) ones, assuming (0, 0) as origin of axes and, incidentally, the position of tree base
#'
#' @param x Abscissa coordinate
#' @param y Ordinate coordinate
#' @return A vector holding angle in degrees and distance in the same unit as \eqn{x} and \eqn{y}
#' @export
toPolar <- function(x, y) {
  d <- sqrt(x^2 + y^2)
  a <- (atan2(x, y) * 180 / pi) %% 360
  c(as.integer(a), d)
}
