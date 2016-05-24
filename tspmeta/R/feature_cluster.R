#' Cluster features.
#'
#' Determines the number of clusters and the mean distances
#' from all cities in a cluster to its centroid.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @param epsilon [\code{numeric(1)}]\cr
#'   Probability in [0,1]. Used to compute the reachability distance
#'   for the underlying \code{\link[fpc]{dbscan}} clustering algorithm.
#' @return [\code{list}].
#' @export
feature_cluster = function(x, epsilon) {
  coords = x$coords
  d = as.vector(x$dists)
  # FIXME: Really strip 0 distances?
	d = d[d > 0]
	q = quantile(d, epsilon)
	# do the clustering
	fit = dbscan(coords, q, showplot = FALSE)
	# skip singleton clusters
	real_clusters = which(fit$cluster > 0)

	cm = fit$cluster[real_clusters]
	coords = coords[real_clusters, , drop = FALSE]
	if (length(cm) > 0) {
		distances = sapply(unique(cm), function(cluster) {
			cluster_coords = coords[cm == cluster, , drop = FALSE]
			centroid = colMeans(cluster_coords)
			apply(cluster_coords, 1, function(point) {
				v = point - centroid
				l2_norm(v)
			})
		})
		distances = unlist(distances)
		res = list(
			number_of_clusters = length(unique(cm)),
			mean_distance_to_centroid = mean(distances)
		)
	} else {
		res = list(
			number_of_clusters = 0,
			mean_distance_to_centroid = NA
		)
	}
	prefix = sprintf("cluster_%02ipct", floor(epsilon * 100))
	names(res) = paste(prefix, names(res), sep = "_")
	res
}