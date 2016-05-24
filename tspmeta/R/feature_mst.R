#' MST features.
#'
#' Construct minimun spanning tree, then calculate
#' the statistics of
#' a) the distances in the MST,
#' b) the depths of all nodes in the MST.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}].
#' @export
feature_mst = function(x) {
	d = x$dists
	# compute spanning tree
	span_tree = spantree(d)
	# depths of MST
	span_tree_depth = spandepth(span_tree)
	# distances within MST
	span_tree_dists = span_tree$dist

	res = c(numvec_feature_statistics(span_tree_depth, "mst_depth"),
					numvec_feature_statistics(span_tree_dists, "mst_dists"))
	res$mst_dists_sum = sum(span_tree_dists) / sum(d)
	res
}
