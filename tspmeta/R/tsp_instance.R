#FIXME: do something when 0 occurs in distances
# Keep in mind, that it is not enough to remove cities!

#' Generates a TSP instance S3 object either from city coordinates.
#'
#' @param coords [\code{matrix}]\cr
#'   Numeric matrix of city coordinates,
#'   rows denote cities.
#' @param dists [\code{dist}]\cr
#'   Optional distance matrix containing the inter-city distances.
#'   If not provided, the (euclidean) distances are computed from the coordinates.
#' @return [\code{\link{tsp_instance}}].
#' @export
tsp_instance = function(coords, dists) {
	assertMatrix(coords, mode = "numeric")
	# FIXME: we need a better check in bbmisc!
	if (!(is.numeric(coords) && nrow(coords) > 3 && ncol(coords) == 2)) {
		stopf("'coords' must be a numeric matrix with at least 4 rows and is currently resctriced to coordinates in 2D! Class = %s, dim = %s",
		  class(coords)[1], collapse(as.character(dim(coords))))
  }
  if (missing(dists)) {
    dists = dist(coords)
  } else {
		checkArg(dists, c("matrix", "dist"))
	}
	dists = as.matrix(dists)
	# FIXME: ugly code
	dists2 = dists
	diag(dists2) = Inf
	if (any(dists2 == 0)) {
		warning("Zero entry detected in TSP dist matrix!")
	}

	structure(list(
		coords = coords,
		dists = dists
	),
  #FIXME: do we limit ourself to euclidean TSP?
  class = c("tsp_instance_euclidean_coords",
            "tsp_instance_symmetric",
            "tsp_instance"))
}

#' Get number of cities in tsp instance.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{integer(1)}].
number_of_cities = function(x) {
	nrow(x$coords)
}

#' Get instance dimensionality (space where coords live).
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{integer(1)}].
instance_dim = function(x) {
	ncol(x$coords)
}

#' Print TSP instance
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param ... [any]\cr
#'   Not used.
print.tsp_instance = function(x, ...) {
  catf("Type                : %c", x$type)
	catf("Coord dimension     : %i", instance_dim(x))
	catf("Number of cities    : %i", number_of_cities(x))
}

#' Plot TSP instance.
#'
#' @param object [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param opt_tour [\code{\link[TSP]{TOUR}}]\cr
#'   TOUR object from package TSP, containing order of cities, tour length and
#'   method name that generated this solution.
#' @param ... [any]\cr
#'   Not used.
#' @return [\code{\link[ggplot2]{ggplot}}].
#' @export
autoplot.tsp_instance = function(object, opt_tour, ...) {
  # draw cities
  coords = data.frame(object$coords)
  names(coords) = c("x", "y")
  pl = ggplot(data = coords, aes_string(x = "x", y = "y"))
  # draw optimal tour
  if (!missing(opt_tour)) {
    # extract order of cities
    opt_tour = as.numeric(opt_tour)
    opt_tour = c(opt_tour, opt_tour[1])
    opt_tour_coords = coords[opt_tour, ]
    pl = pl + geom_path(data = opt_tour_coords, colour = "tomato")
  }
  pl = pl + geom_point(colour = "tomato")
  #pl = pl + xlim(c(0,1))
  #pl = pl + ylim(c(0,1))
  pl
}


# FIXME: as. with dot?
#' Convert to TSP instance object of package TSP.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{\link[TSP]{TSP}}].
#' @export
as_TSP = function(x) {
  TSP(x$dists)
}
