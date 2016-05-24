causal.effect <-
function(y, x, z = NULL, G, expr = TRUE, simp = TRUE) {
  if (!is.dag(observed.graph(G))) stop("Graph 'G' is not a DAG")
  to <- topological.sort(observed.graph(G))
  to <- get.vertex.attribute(G, "name")[to]
  G.unobs <- unobserved.graph(G)
  G.Adj <- as.matrix(get.adjacency(G.unobs))
  if (length(setdiff(y, to)) > 0) stop("Set 'y' contains variables not present in the graph.")
  if (length(setdiff(x, to)) > 0) stop("Set 'x' contains variables not present in the graph.")
  if (length(z) > 0) {
    if (length(setdiff(z, to)) > 0) stop("Set 'z' contains variables not present in the graph.")
  }
  res <- probability()
  if (length(z) == 0) { 
    res <- id(y, x, probability(), G, to)
    # res <- organize.terms(res)
  } else { 
    res <- idc(y, x, z, probability(), G, to)
    # res <- organize.terms(res)
    res2 <- res
    res2$sumset <- union(res2$sumset, y)
    res$fraction <- TRUE
    res$divisor <- res2
  }
  if (simp) {
    to.u <- topological.sort(G.unobs)
    to.u <- get.vertex.attribute(G.unobs, "name")[to.u]
    res <- deconstruct(res, probability(recursive = TRUE, children = list(), fraction = TRUE, divisor = probability(recursive = TRUE, children = list())))
    res <- parse.expression(res, to, G.Adj)
    if (length(res$divisor$children) > 0) {
      res$divisor <- parse.expression(res$divisor, to, G.Adj)
      if (zero.children(res$divisor)) {
        res$fraction <- FALSE
        res$divisor <- list()
      }
    } else {
      res$fraction <- FALSE
      res$divisor <- list()
    }
    res <- deconstruct(res, probability(recursive = TRUE, children = list(), fraction = TRUE, divisor = probability(recursive = TRUE, children = list())))
    res <- parse.deconstruct(res)
    res <- simplify.expression(res, G.unobs, to.u)
  }
  if (expr) res <- get.expression(res)
  return(res)
}
