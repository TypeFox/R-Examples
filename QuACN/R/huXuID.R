huXuID <- function(g, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  if (is.null(deg))
    deg <- graph::degree(g)

  n <- nodes(g)
  nd <- .nodeDataVector(g, "atom->number")
  e <- edges(g)
  ed <- .edgeDataMatrix(g, "bond")

  deg[deg == 0] <- 0.5
  Z <- sapply(n, function(v) as.integer(nd[[v]]))
  hxdeg <- deg * sqrt(Z)

  weightfunc <- function(i, from, to) {
    sqrt(
     (as.numeric(ed[[from, to]]) / (i + 1)) *
     (1 / (hxdeg[[from]] * hxdeg[[to]]))
    )
  }

  AID <- sapply(n, function(v) .weightedPathSum(e, v, weightfunc, unit=0, rmdupls=FALSE))

  sum(AID ^ 2)
}
