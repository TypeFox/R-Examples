hasCycles <- function(A, verbose = FALSE){ 
  G <- graph.adjacency(A,mode="directed",weighted="a")
  cyc <- !is.dag(G)
  if(verbose) if(cyc) cat('Graph is cyclic\n') else cat('Graph is acyclic\n')
  cyc
}