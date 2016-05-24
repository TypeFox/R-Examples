# Select cluster centers to form depots.
#
# @param n.depots [integer(1)]
#   Number of depots. Currently at most 2.
# @param cluster.centers [\matrix]
#   Matrix containing the coordinates of the cluster centers.
# @param distances [list]
#   See return value of computeDistancesToNearestClusterCenters.
# @return [matrix]
#   Coordinates of the depots.
buildDepots = function(n.depots, cluster.centers, distances) {
  n.cluster = nrow(cluster.centers)
  # get first depot randomly
  depot.1.idx = sample(seq(n.cluster), 1L)
  depot.coordinates = cluster.centers[depot.1.idx, , drop = FALSE]
  if (n.depots == 2L) {
    depot.2.idx = distances$max.distance.idx[depot.1.idx]
    depot.coordinates = rbind(depot.coordinates, cluster.centers[depot.2.idx, , drop = FALSE])
  }
  return(depot.coordinates)
}
