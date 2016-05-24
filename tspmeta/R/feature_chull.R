#' Convex hull features.
#'
#' Determines the area of the convex hull and the
#' ratio of the cities which lie on the convex hull in
#' the euklidean space.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}].
#' @export
feature_chull = function(x) {
    coords = x$coords
    hull = chull(coords[, 1], coords[, 2])
    area = areapl(coords[hull, ])
    list(
        chull_area = area,
	   chull_points_on_hull = length(hull) / nrow(coords)
    )
}
