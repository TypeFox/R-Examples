informationBondIndex <- function(g) {
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  m <- numEdges(g)

  bonds <- .edgeDataMatrix(g, "bond")
  bonds <- bonds[upper.tri(bonds)]
  eqcls <- as.numeric(table(bonds[bonds != 0]))

  m * log2(m) - sum(eqcls * log2(eqcls))
}
