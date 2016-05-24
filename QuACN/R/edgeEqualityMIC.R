edgeEqualityMIC <- function(g, deg=NULL) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  if (is.null(deg))
    deg <- graph::degree(g)

  conns <- .edgeApply(g, function(from, to) deg[[from]] * deg[[to]],
                      dupls=FALSE)
  cls <- table(conns)
  p <- as.numeric(cls / numEdges(g))

  -sum(p * log2(p))
}
