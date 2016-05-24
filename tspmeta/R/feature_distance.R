#' Distance features.
#'
#' Computes different statistics describing the distribution of
#' pairwise distances between cities.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}]
#'   List of statistics describing the distribution of distances.
#' @export
feature_distance = function(x) {
	size = number_of_cities(x)
  dd = x$dists
  freq = table(as.vector(dd))
  mode_freq = as.numeric(max(freq)) # as.numeric to strip name attribute
  mode_values = as.numeric(names(freq[freq == mode_freq]))
  mode_quantity = length(mode_values) / length(dd)
  res = list(
       distance_distances_shorter_mean_distance = sum(dd < mean(dd)) / length(dd),
       distance_distinct_distances = length(unique(dd)) / length(dd),
       distance_mode_frequency = mode_freq,
       distance_mode_quantity = mode_quantity,
       distance_mode_mean = mean(mode_values),
       distance_mean_tour_length = 2 / (size - 1) * sum(dd),
       distance_sum_of_lowest_edge_values = sum(sort(as.dist(dd), partial = size)[1:size])
       )
  c(res, numvec_feature_statistics(dd, "distance"))
}