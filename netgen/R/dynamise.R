#' Add dynamic arrival times to nodes.
#'
#' Some variants of the Vehicle Routing Problem (VRP) consider static as well
#' as dynamic customers (nodes). This function takes a \code{Network} and
#' dynamises it, i. e., it adds dynamic arrival times to the customers via a
#' Poisson process.
#'
#' @template arg_network
#' @param n.dynamic [integer(1) | NULL]
#'   Number of nodes, which should become dynamic. Ignored if \code{dyn.customers.ratio}
#'   is not \code{NULL}.
#' @param dyn.customers.ratio [numeric(1) | NULL]
#'   Ratio of dynamic nodes. If this is set to a numeric value in (0, 1), the
#'   parameter \code{n.dynamic} is ignored.
#' @param arrival.limit [numeric(1)]\cr
#'   Maximal arrival time.
#' @return [\code{Network}]
#'   Modified network (now has an additional list element 'arrival.times') and the
#'   ratio of dynamic customers as an attribute.
#' @examples
#' x = generateClusteredNetwork(n.points = 100L, n.cluster = 4L, upper = 100, n.depots = 2L)
#' x = dynamise(x, dyn.customers.ratio = 0.3, arrival.limit = 400)
#' print(x)
#' @seealso \code{\link{generateRandomNetwork}}, \code{\link{generateClusteredNetwork}},
#' \code{\link{generateGridNetwork}}
#' @export
dynamise = function(x, n.dynamic = NULL, dyn.customers.ratio = NULL, arrival.limit) {
  assertClass(x, "Network")
  assertNumber(arrival.limit, lower = 1, finite = TRUE)

  if (is.null(n.dynamic) && is.null(dyn.customers.ratio)) {
    stopf("Either n.dynamic or dyn.customers.ratio must be set.")
  }

  if (!is.null(n.dynamic)) {
    assertInteger(n.dynamic, lower = 1, upper = getNumberOfNodes(x))
  }

  n.customers = getNumberOfNodes(x)
  # ignore n.dynamic if ratio is set
  if (!is.null(dyn.customers.ratio)) {
    assertNumber(dyn.customers.ratio, lower = 0, upper = 1)
    n.dynamic = ceiling(n.customers * dyn.customers.ratio)
  }
  # determine rate of arrival times
  rate = arrival.limit / n.dynamic

  # add preliminary request times
  x$arrival.times = rep(0, n.customers)

  # randomly select n.dynamic customers
  idx.dyn = sample(1:n.customers, replace = FALSE, size = n.dynamic)

  # sample arrival times according to Poisson process
  x$arrival.times[idx.dyn] = sampleArrivalTimes(n.dynamic, arrival.limit, 1 / rate)

  attr(x, "dyn.customers.ratio") = dyn.customers.ratio
  return(x)
}

# Helper function which sample the arrival times.
#
# Internally the function generates a sequence of 'request times' and builds
# the cumulative sum to get arrival times. It potentially repeats the process
# multiple times to ensure, that the last arrival time is lower than the
# arrival.limit.
#
# @param n [integer(1)]
#   Number of arrival times.
# @param arrival.limit [numeric(1)]
#   Limit for the arrival times.
# @param rate [numeric(1)]
#   Rate of the Poisson process, i. e., the expected value of the
#   exponential distribution.
# @return [numeric(n)]
#   Arrival times.
sampleArrivalTimes = function(n, arrival.limit, rate) {
  arrival.times = cumsum(rexp(n, rate = rate))
  while (any(arrival.times > arrival.limit)) {
    arrival.times = cumsum(rexp(n, rate = rate))
  }
  arrival.times
}
