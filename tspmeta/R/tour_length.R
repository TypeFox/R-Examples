#' @export
tour_length.tsp_instance = function(x, order, ...) {
  tour_length(as_TSP(x), ...)
}
