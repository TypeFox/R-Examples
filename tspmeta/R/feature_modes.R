#' Modes of edge cost distribution feature.
#'
#' Includes the number of modes of the edge cost distribution.
#'
#' @param x [\code{\link{tsp_instance}}]\cr
#'   TSP instance.
#' @return [\code{list}]
#'   List containing (estimated) number of modes.
#' @export
feature_modes = function(x) {
  dists = as.numeric(x$dists)

  intdens = function(a, b) {
    mean(y[a:b]) * (d$x[b] - d$x[a])
  }
  
  d = density(dists)
  y = d$y
  n = length(y)
  minidx = c(1, which(y[2:(n-1)] < y[1:(n-2)] & y[2:(n-1)] < y[3:n]), n + 1)
  modemass = sapply(1:(length(minidx) - 1),
                     function(i) intdens(minidx[i], minidx[i + 1] - 1))
  list(modes_number = sum(modemass > 0.01))
}
