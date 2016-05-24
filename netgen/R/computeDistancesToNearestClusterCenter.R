# Helper function for the computation of the distances to the nearest/farthest
# neighbour cluster center.
#
# @param cluster.centers [matrix]
#   Numeric matrix. Each row corresponds to the coordinates of a cluster center.
# @return [list]
#   List of two elements.
#   \describe{
#     \item{min.distance}{Numeric vector of the distance to the nearest neighbor}
#     \item{min.distance.idx}{Integer vector. The entry is the row id of the nearest
#     neighbor cluster center coordinates.}
#     \item{max.distance}{Numeric vector of the distance to the farthest neighbor}
#     \item{max.distance.idx}{Integer vector. The entry is the row id of the farthest
#     neighbor cluster center coordinates.}
#   }
computeDistancesToNearestClusterCenter = function(cluster.centers) {
  n = nrow(cluster.centers)
  min.distance.idx = numeric(n)
  max.distance.idx = numeric(n)
  min.distance = numeric(n)
  max.distance = numeric(n)
  for (i in seq(n)) {
    # compute pairwise euclidean distances
    distances = apply(cluster.centers, 1, function(x) {
      sqrt(sum((x - cluster.centers[i, ])^2))
    })
    # since the distance to itself is always zero we set this to Infinity here
    distances[i] = Inf
    min.distance.idx[i] = which.min(distances)
    min.distance[i] = min(distances)
    distances[i] = 0
    max.distance[i] = max(distances)
    max.distance.idx[i] = which.max(distances)
  }
  return(list(
    min.distance = min.distance,
    min.distance.idx = min.distance.idx,
    max.distance = max.distance,
    max.distance.idx = max.distance.idx
  ))
}
