edgeDeletedSubgraphs <- function(gs) {
  if (class(gs) != "list")
    gs <- list(gs)

  raw_result <- lapply(gs, function(g) {
    if (class(g)[1] == "graphNEL")
      M <- adjacencyMatrix(g)
    else if (class(g)[1] == "matrix")
      M <- g
    else
      stop("'gs' must be a list of 'graphNEL' objects or adjacency matrices")

    E <- which(upper.tri(M) & (M == 1), TRUE)
    mapply(function(from, to) {
      M_copy <- M
      M_copy[from, to] <- 0
      M_copy[to, from] <- 0
      M_copy
    }, E[, "row"], E[, "col"], SIMPLIFY=FALSE)
  })

  unique(unlist(raw_result, FALSE))
}
