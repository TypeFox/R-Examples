## finds ancestors of some given node(s)

findAncestors <- function(G, node, graph = FALSE) {
  n <- nrow(G) 
  dag <- matrix(0, n, n )
  dag[ G == 2 & t(G) == 3 ] <- 1
  isAnc <- transitiveClosure(dag)
  
  anc <- as.vector( isAnc[, node] )
  anc <- which( anc > 0 )
  Ganc <- G[c(node, anc), c(node, anc)]

  if (graph == TRUE) {
    if ( length(anc) > 0 ) {
      if ( requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE ) {
        Ganc[ Ganc != 2 ] <- 0
        g <- as( Ganc, "graphNEL" )
        plot(g, main = paste("Completed partially directed graph with ancestors of ", node ) )
      } else {
        warning('In order to plot the generated network, package Rgraphviz is required.')
      }
    } 
  }
  
  list(isAnc = isAnc, Ganc = Ganc, anc = anc)   
}