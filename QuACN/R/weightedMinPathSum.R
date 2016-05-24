.weightedMinPathSum <- function(g, weightfunc) {
  stopifnot(.validateGraph(g))
  n <- nodes(g)
  ig <- .G2IG(g)

  sum(sapply(1:length(n), function(v) {
    sp <- igraph::get.all.shortest.paths(ig, from=v)$res
    sum(sapply(sp, function(path) {
      if (length(path) == 1)
        1
      else
        prod(sapply(1:(length(path) - 1), function(i) {
          weightfunc(i, n[[path[[i]]]], n[[path[[i + 1]]]])
        })) / 2
    }))
  }))
}
