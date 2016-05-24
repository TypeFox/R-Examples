#' Rescale coords of TSP instance to \eqn{[0,1]^2}.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param coords [\code{matrix}]\cr
#'   Numeric matrix of city coordinates,
#'   rows denote cities.
#' @return [\code{matrix}] for \code{rescale_coords} and \code{tsp_instance} for \code{rescale_instance}.
#'   Numeric matrix of scaled city coordinates.
#' @export
rescale_instance = function(x) {
    assertClass(x, "tsp_instance")
    x$coords = rescale_coords(x$coords)
    x$dists = as.matrix(dist(x$coords))
    return(x)
}

#' @export
#' @rdname rescale_instance
rescale_coords = function(coords) {
    min = apply(coords, 2, min)
    max = apply(coords, 2, max)
    t((t(coords) - min) / (max - min))
}


