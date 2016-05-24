#' Rotate a matrix of 2D coordinates
#'
#' @param coords [\code{matrix}]\cr
#'   Numeric matrix of 2D coordinates to rotate
#' @param angle [\code{numeric(1)}]\cr
#'   Angle by which to rotate the coordinates. In radians.
#' @param center [\code{matrix}]\cr
#'   Center around which to rotate the coordinates.
#'
#' @return A matrix of rotated coordinates.
rotate_coordinates = function(coords, angle, center) {
  R = matrix(c(cos(angle), -sin(angle),
                sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
  t((R %*% (t(coords) - center)) + center)
}

#' Rotate the cities of a TSP instance around a point.
#'
#' @param instance [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param angle [\code{numeric(1)}]\cr
#'   Angle by which to rotate the coordinates. In radians.
#' @param center [\code{numeric}]\cr
#'   Point around which to rotate the cities. If missing, defaults to
#'   the center of mass of the cities.
#'
#' @return [\code{\link{tsp_instance}}] New TSP instance.
#'
#' @export
rotate_instance = function(instance, angle, center) {
  assertClass(instance, "tsp_instance")
  assertNumber(angle, na.ok = FALSE)
  if (missing(center)) {
    center = center_of_mass(instance)
  }
  assertNumeric(center, len = 2L, any.missing = FALSE)
  instance$coords = rotate_coordinates(instance$coords, angle, center)
  instance
}

#' Calculate rotation angle such that the main axis through the
#' cities is aligned with the X axis.
#'
#' @param instance [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#'
#' @return [\code{numeric(1)}]
#'
#' @export
normalization_angle = function(instance) {
  assertClass(instance, "tsp_instance")
  C = cov(instance$coords)
  main_axis = eigen(C)$vectors[, 1]
  angle = -atan2(main_axis[2], main_axis[1])
  angle
}

#' Normalize an instance w.r.t. its rotation.
#'
#' Normalization is performed by aligning the main axis of the cities
#' with the X axis.
#'
#' @param instance [\code{tsp_instance}]\cr
#' @return A rotated \code{tsp_instance}.
#' @seealso \code{\link{normalization_angle}}
#' @export
normalize_rotation = function(instance) {
  ## Rotate by an additional 45deg so that the axis is aligned with
  ## the line with slope 1. This is mainly for cosmetic reasons when
  ## the instances are normalized to lie within the unit square.
  angle = normalization_angle(instance) + pi / 4
  cg = center_of_mass(instance)
  i1 = rotate_instance(instance, angle, cg)
  i2 = rotate_instance(instance, angle + pi, cg)

  ## Decide which instance is 'better' by counting the number of
  ## points to the right of the center of gravity.
  if (sum(i1$coords[, 1] > cg[1]) > sum(i2$coords[, 1] > cg[1]))
    i1
  else
    i2
}
