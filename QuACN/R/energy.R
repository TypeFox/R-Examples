energy <- function(g) {
 
  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  M <- adjacencyMatrix(g)
  EV <- as.double(eigen(M, only.values=TRUE)$values)
  sum(abs(EV))
}
