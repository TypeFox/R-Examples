#####################
##### Plot a graph using the adjacency matrix
#####
####################

plota <- function(G) {
  if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE) {
    if ( sum(lower.tri(G - t(G) ) ) == 0 ) { ## no orientations
      g <- new("graphAM", adjMat = G, edgemode = "undirected")
      plot(g, main = "Association network graph")
    } else {
      G[ G == 1 ] <- 2
      G[G != 2] <- 0
      g <- as( G, "graphNEL" )
      plot(g, main = paste("Completed partially directed graph" ) )
    }
  } else {
    warning('In order to plot the generated network, package Rgraphviz is required.')
  }
}
