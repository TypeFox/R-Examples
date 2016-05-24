# Helper function for distributing desired number of points
# more or less equally to the different clusters.
#
# @param n.cluster [\code{integer(1)}]\cr
#   Desired number of clusters.
# @param n.points [\code{integer(1)}]\cr
#   Number of points for the instance.
# @param strategy [\code{character(1)}]\cr
#   Strategy used to determine number of points per cluster.
#   See \code{getPointDistributionStrategy} for a list of available strategies.
# @return [\code{integer}]
#   Vector of length \code{length(n.cluster)} containing the number of
#   points assigned to this cluster.
determineNumberOfPointsPerCluster = function(n.cluster, n.points, strategy = "equally.distributed") {
  if (strategy == "equally.distributed") {
    getEquallyDistributedIntegerPartition(n.points, n.cluster)
  } else if (strategy == "random.partition") {
    getRandomIntegerPartition(n.points, n.cluster)
  }
}

# Computes integer partition.
#
# Integer n should be partitioned in k integers, i1, ...,ik that way, that
# the i1 + ... + ik = n and no ij is zero and i1 approx i2 approx ... approx ik.
#
# @param n [integer(1)]
#   Integer we want to get an integer partiton for.
# @param k [integer(1)]
#   Number of partitions.
# @return [integer(k)]
getEquallyDistributedIntegerPartition = function(n, k) {
  # distribute equally over the clusters
  m = floor(n / k)
  partition = rep(m, k)
  # n * n.cluster might be lower than n.points. Add the remaining points to
  # a randomly chosen cluster
  # FIXME: we might want to implement different strategies here
  m.diff = n - m * k
  # sample random cluster
  idx = sample(1:k, 1)
  partition[idx] = m + m.diff
  return(partition)
}

# Computes integer partition.
#
# Integer n should be partitioned in k integers, i1, ...,ik that way, that
# the i1 + ... + ik = n and no ij is zero.
#
# @param n [integer(1)]
#   Integer we want to get an integer partiton for.
# @param k [integer(1)]
#   Number of partitions.
# @return [integer(k)]
#FIXME: this does not work very well. First numbers are big, later small,
#occasionally negative values in the end.
#FIXME: this is currently disabled, because it is buggy! Next release.
getRandomIntegerPartition = function(n, k) {
  s = n
  rest = s
  # set up storage for partition
  partition = integer(k)
  for (i in seq(k)) {
    # if the last partition is handled, we need to assign the 'rest' to it
    if (i == k) {
      partition[i] = rest
      return(partition)
    }
    # otherwise sample a random number of maximal size
    partition[i] = sample(1:(floor(rest / 2)), size = 1L)
    rest = rest - partition[i]
  }
  return(partition)
}

#' Returns the available strategies for distributing points around clusters.
#'
#' @return [\code{character}]
#' @export
#FIXME: we do not export "random.partition", since there is a bug in the
# implementation so far.
getPointDistributionStrategies = function() {
  return(c("equally.distributed")) #, "random.partition"))
}
