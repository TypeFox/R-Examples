.edgeDataMatrix <- function(g, att) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  if (is.null(edgeDataDefaults(g)[[att]]))
    stop(paste("edge attribute", att, "not set"))

  raw <- edgeData(g)
  n <- nodes(g)
  m <- matrix(0, nrow=length(n), ncol=length(n), dimnames=list(n, n))
  for (e in names(raw)) {
    vert <- strsplit(e, "\\|")[[1]]
    m[[vert[[1]], vert[[2]]]] <- raw[[e]][[att]]
  }
  m
}
