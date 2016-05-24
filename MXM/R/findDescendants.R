## finds descenadants of some given node(s)

findDescendants <- function(G, node, graph = FALSE) {
  n <- nrow(G) 
  dag <- matrix(0, n, n)
  dag[ G == 3 & t(G) == 2 ] <- 1
  isDesc <- transitiveClosure(dag)
  
  desc <- as.vector( isDesc[, node] )
  desc <- which( desc > 0 )
  Gdesc <- G[c(node, desc), c(node, desc)]
  
  if (graph == TRUE) {
    if ( length(desc) > 0 ) {
      if ( requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE ) {
        Gdesc[ Gdesc != 2 ] <- 0
        g <- as( Gdesc, "graphNEL" )
        plot(g, main = paste("Completed partially directed graph with ancestors of ", node ) )
      } else {
        warning('In order to plot the generated network, package Rgraphviz is required.')
      }
    } 
  }
  
  list(isDesc = isDesc, Gdesc = Gdesc, desc = desc)     
}