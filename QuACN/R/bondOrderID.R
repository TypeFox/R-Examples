bondOrderID <- function(g) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))

  e <- edges(g)
  ed <- .edgeDataMatrix(g, "bond")
  .weightedPathSum(e, NULL, function(i, from, to) as.numeric(ed[[from, to]]))
}
