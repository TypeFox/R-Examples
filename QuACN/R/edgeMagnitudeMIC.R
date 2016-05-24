edgeMagnitudeMIC <- function(g, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  if (is.null(deg))
    deg <- graph::degree(g)

  conns <- .edgeApply(g, function(from, to) 1/sqrt(deg[[from]] * deg[[to]]),
                      dupls=FALSE)
  randic <- sum(conns)
  p <- conns / randic

  -sum(p * log2(p))
}
