aux.effect <-
function(y, x, z, G, expr = TRUE, simp = TRUE) {
  if (!is.dag(observed.graph(G))) stop("Graph 'G' is not a DAG")
  to <- topological.sort(observed.graph(G))
  to <- get.vertex.attribute(G, "name")[to]
  if (length(setdiff(y, to)) > 0) stop("Set 'y' contains variables not present in the graph.")
  if (length(setdiff(x, to)) > 0) stop("Set 'x' contains variables not present in the graph.")
  if (length(z) > 0) {
    if (length(setdiff(z, to)) > 0) stop("Set 'z' contains variables not present in the graph.")
  }
  if (length(z) == 0) { 
    cat("Reverting to ordinary identifiability \n")
    res <- id(y, x, probability(), G, to)
  } else { 
    res <- zid(y, x, z, NULL, NULL, probability(), G, to)
  }
  if (simp) {
    G.unobs <- unobserved.graph(G)
    to.u <- topological.sort(G.unobs)
    to.u <- get.vertex.attribute(G.unobs, "name")[to.u]
    res <- simplify.expression(res, G.unobs, to.u)
  }
  res <- organize.terms(res)
  if (expr) res <- get.expression(res)
  return(res)
}
