leaderCluster = function(points, radius, weights = rep(1, nrow(points)), max_iter = 10L,
                          distance = c("Lp", "L1", "L2", "Linf", "haversine"), p = 2)
{
  type_code = match(distance[[1]], c("Lp", "L1", "L2", "Linf", "haversine")) - 1L
  if (!is.matrix(points)) points = as.matrix(points)
  if (is.na(type_code))
  {
    stop(paste0("The distance method '", distance[[1]], "' is not supported."))
  }
  if (p <= 0) stop(paste0("The input p='", p, "' is invalid. Must be > 0."))
  if (type_code == 0L & p == 1) type_code = 1L
  if (type_code == 0L & p == 2) type_code = 2L
  if (type_code == 0L & p == Inf) type_code = 3L
  if (!is.finite(p)) p = 2
  if (type_code == 4L & ncol(points) != 2L)
  {
    stop("Haversine distance requires two dimensional input.")
  }

  cluster_id = rep(0, nrow(points))

  iter = 1L
  index = 1L:nrow(points)
  ord = order(weights, decreasing=TRUE)
  points = points[ord,,drop=FALSE]
  index = index[ord]

  while (iter <= max_iter )
  {

    out = .C("leader_cluster",
              delta = as.double(radius),
              points = as.double(points),
              weights = as.double(weights),
              cluster_id = as.integer(cluster_id),
              nrow = nrow(points),
              ncol = ncol(points),
              type = type_code,
              p = p,
              PACKAGE = "leaderCluster")$cluster_id + 1

    if(all(out == cluster_id)) break

    cluster_id = out
    ord = order(out)
    points = points[ord,,drop=FALSE]
    index = index[ord]
    iter = iter + 1
  }

  num_clusters = max(cluster_id)
  cluster_centroids = matrix(NA, ncol=ncol(points), nrow=num_clusters)
  for (i in 1:ncol(points)) cluster_centroids[,i] = tapply(points[,i], cluster_id, mean)

  return(list(cluster_id = out[match(1:length(index), index)],
              cluster_centroids = cluster_centroids,
              num_clusters = num_clusters,
              iter = iter))
}
