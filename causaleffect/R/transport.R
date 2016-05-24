transport <-
function(y, x, z = NULL, D, expr = TRUE, simp = TRUE) {
  if (!is.dag(observed.graph(D))) stop("Selection diagram 'D' is not a DAG")
  to <- topological.sort(observed.graph(D))
  to <- get.vertex.attribute(D, "name")[to]
  if (length(setdiff(y, to)) > 0) stop("Set 'y' contains variables not present in the graph.")
  if (length(setdiff(x, to)) > 0) stop("Set 'x' contains variables not present in the graph.")
  if (length(z) > 0) {
    if (length(setdiff(z, to)) > 0) stop("Set 'z' contains variables not present in the graph.")
  }
  res <- probability()
  if (length(z) == 0) res <- sid(y, x, probability(star = TRUE), D, to)
  else {
    res <- trz(y, x, probability(star = TRUE), z, NULL, D, to)
  }
  if (simp) {
    D.unobs <- unobserved.graph(D)
    to.u <- topological.sort(D.unobs)
    to.u <- get.vertex.attribute(D.unobs, "name")[to.u]
    res <- simplify.expression(res, D.unobs, to.u)
  }
  res <- organize.terms(res)
  if (expr) res <- get.expression(res)
  return(res)
}
