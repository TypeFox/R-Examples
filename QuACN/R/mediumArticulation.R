mediumArticulation <- function(g) {

  if (class(g)[1] != "graphNEL")
    stop("'g' must be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  n <- numNodes(g)
  m <- numEdges(g)
  R_Clique <- 2*log(n-1)
  R_path <- 2*((n-2)/(n-1))*log(2)
  I_Clique <- log(n/(n-1))
  I_path <- log(n-1) - (n-3)/(n-1)*log(2)
  deg <- graph::degree(g)
  M <- adjacencyMatrix(g)

  edges <- which(upper.tri(M, TRUE) & (M == 1), TRUE)
  from <- edges[, "row"]
  to <- edges[, "col"]

  R <- sum(mapply(function(j, k) log(deg[j] * deg[k]), from, to)) / m
  I <- sum(mapply(function(j, k) log((2 * m) / (deg[j] * deg[k])), from, to)) / m

  MA_R <- (4*(R-R_path)*(R_Clique-R)) / ((R_Clique-R_path)^2)
  MA_I <- (4*(I_path-I)*(I-I_Clique)) / ((I_path-I_Clique)^2)

  MA_R * MA_I
}
