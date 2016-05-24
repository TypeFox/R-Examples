minConnectivityID <- function(g, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(deg))
    deg <- graph::degree(g)

  .weightedMinPathSum(g, weightfunc=function(i, from, to) (1 / sqrt(deg[[from]] * deg[[to]])))
}
