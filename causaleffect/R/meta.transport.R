meta.transport <-
function(y, x, D, expr = TRUE, simp = TRUE) {
  d <- length(D)
  to <- lapply(D, function(x) topological.sort(observed.graph(x)))
  to <- lapply(1:d, function(x) get.vertex.attribute(D[[x]], "name")[to[[x]]])
  for (i in 1:d) {
    if (!is.dag(observed.graph(D[[i]]))) stop("Selection diagram 'D[", i, "]' is not a DAG")
    if (length(setdiff(y, to[[i]])) > 0) stop("Set 'y' contains variables not present in the diagrams.")
    if (length(setdiff(x, to[[i]])) > 0) stop("Set 'x' contains variables not present in the diagrams.")
  }
  res <- usid(y, x, probability(star = TRUE), D, to)
  if (simp) {
    D.unobs <- lapply(D, function(x) unobserved.graph(x))
    to.u <- lapply(D.unobs, function(x) topological.sort(x))
    to.u <- lapply(1:d, function(x) get.vertex.attribute(D.unobs[[x]], "name")[to.u[[x]]])
    res <- simplify.meta.expression(res, D.unobs, to.u)
  }
  res <- organize.terms(res)
  if (expr) res <- get.expression(res)
  return(res)
}
