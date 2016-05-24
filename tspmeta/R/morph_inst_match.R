# Component-by-component convex combination of the coordinates of 2 TSP instances.
#
# Combination is alpha * x + (1 - alpha) * y.
#
# @param x [\code{matrix}]\cr
#   Numeric matrix of city coordinates.
# @param y [\code{matrix}]\cr
#   Numeric matrix of city coordinates.
# @param alpha [\code{numeric(1)}]\cr
#   Coefficient alpha for convex combination
# @return [\code{matrix}]
#   New coordinates.
convex_combination = function(x_coords, y_coords, alpha) {
  alpha * x_coords + (1 - alpha) * y_coords
}

#' Morphing (convex-combination) of two instances with parameter alpha.
#'
#' Pairs of cities are matched in a greedy fashion, see \code{\link{greedy_point_matching}}.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#"   TSP instance.
#' @param y [\code{\link{tsp_instance}}]\cr
#"   TSP instance.
#' @param alpha [\code{numeric(1)}]\cr
#'   Coefficient alpha for convex combination.
#' @return [\code{\link{tsp_instance}}]
#'   Morphed TSP instance.
#' @export
#' @examples
#' x = random_instance(10)
#' y = random_instance(10)
#' z = morph_instances(x, y, 0.5)
#' autoplot(x)
#' autoplot(y)
#' autoplot(z)
morph_instances = function(x, y, alpha) {
    assertClass(x, "tsp_instance")
    assertClass(y, "tsp_instance")
    assertNumber(alpha, lower = 0, upper = 1, na.ok = FALSE)

    points_m = greedy_point_matching(x, y)
    coords = convex_combination(x$coords, y$coords[points_m[, 2]], alpha)

    tsp_instance(coords = coords)
}
