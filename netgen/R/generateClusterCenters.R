# Helper function to generate cluster centres.
#
# The function generates n cluster centers by generating a LHS design of size
# n in the unit-cube.
#
# @param n.cluster [\code{integer(1)}]\cr
#   Number of clusters.
# @param n.dims [\code{integer(1)}]\cr
#   Number of dimensions. Default is 2.
# @param lower [\code{numeric(1)}]\cr
#   Lower box constaint for cube.
# @param upper [\code{numeric(1)}]\cr
#   Upper box constaint for cube.
# @return [\code{matrix}]
#   Cluster center matrix. Each row contains the coordinates of one cluster center.
generateClusterCenters = function(
  n.cluster = 5L, n.dims = 2L,
  generator = lhs:::maximinLHS,
  lower = 0, upper = 1) {
  # FIXME: introduce use of min.dist.to.bounds
  cc = generator(n.cluster, n.dims)

  # "stretch design"
  cc = lower + (upper - lower) * cc
  return(cc)
}
