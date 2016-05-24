############################
############################
##### Return and plot if requested the neighbours of an undirected graph
#####
############################
############################

nei <- function(G, node, graph = TRUE) {
  ## G is the adjacency matrix of an UN-DIRECTED graph
  ## node is either one number between 1 and the number of nodes
  ## or more numbers, corresponding to two or more nodes
  ## it is a node whose neighbors you want to find
  ## it can also be a vector with more than one nodes
  
  if ( is.null( colnames(G) ) ) {
    p <- ncol(G)
    nama <- paste("X", 1:p, sep = "")
  } else  nama <- colnames(G)
   
  n <- length(node)
  if ( n == 1 ) {
    ind = which( G[node, ] == 1 ) 
    if ( length(ind) == 0 ) {
      graph <- FALSE
      geit <- paste( "The chosen node has no neighbours" );
    } else {
      ind <- as.vector( ind[ ind>0 ] )
      geit <- c(node, ind)
      names(geit) <- nama[ c(node, ind) ]
      Gnei <- G[ geit, geit ]
    }

  } else if ( n == 2 ) { 
    geit <- undir.path(G, node[1], node[2])
    geit <- unique( geit[geit > 0] )
    names(geit) <- nama[geit]
    Gnei <- G[geit, geit]    

  } else if ( n > 2 ) {
    posa <- t( combn(node, 2) )
    geit <- list()
    for ( i in 1:nrow(posa) ) {
      geiton <- undir.path(G, posa[i, 1], posa[i, 2])
      geit[[ i ]] <- unique( geiton[geiton > 0] )
    }
    geit <- unlist(geit)   
    geit <- unique( geit[geit > 0] )
    Gnei <- G[geit, geit]
    names(geit) <- nama[geit]
  }

    if (graph == TRUE) {
      if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE) {
        am.graph <- new("graphAM", adjMat = Gnei, edgemode = "undirected")
        plot(am.graph, main = "Subgraph of association network")
      } else {
        warning('In order to plot the generated network, package Rgraphviz is required.')
      }
    }
  geit
}