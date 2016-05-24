.edgeApply <- function(g, f, dupls=TRUE, simplify=TRUE) {
  if (class(g)[1] == "graphNEL")
    M <- adjacencyMatrix(g)
  else
    M <- g

  cond <- (M == 1)
  if (!dupls)
    cond <- cond & upper.tri(M, TRUE)

  edges <- which(cond, TRUE)
  mapply(f, edges[, "row"], edges[, "col"], SIMPLIFY=simplify)
}
