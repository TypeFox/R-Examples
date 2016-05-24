############################
############################
##### Return and plot if requested the neighbours of an undirected graph
#####
############################
############################

mb <- function(G, node, graph = FALSE) {
  ## G is the adjacency matrix of an UN-DIRECTED graph
  ## node is a number between 1 and the number of nodes
  ## it is a node whose neighbors you want to find
  ## it can also be a vector with more than one nodes
  
  if ( is.null( colnames(G) ) ) {
    p <- ncol(G)
    nama <- paste("X", 1:p, sep = "")
  } else  nama <- colnames(G)
  
  parents <- which( G[node, ] == 1 & t(G[node, ] != 1) )
  relatives <- which( ( G[node, ] == 1 & G[, node] == 1 ) )
  children <- which( G[node, ] == 2 )
  le <- length(children) 
  spouses <- list()
  if (le > 0) {
    for (i in le) {
      spousa <- which( G[children, ] == 1 )
      spouses[[ i ]] <- setdiff(spousa, node) 
    }
  }
  
  names(spouses) <- nama[ children ]
  spo <- unlist(spouses)
  blanket <- c(parents, children, spo, relatives)
  blanket <- unique(blanket)
  
  if ( length(relatives) == 0 ) {
    graph <- FALSE
  } else {  
    Grel <- G[c(node, relatives), c(node, relatives)]
    aa <- which(Grel == 1 & t(Grel) == 1, arr.ind = TRUE)
    Grel[aa] <- 2
    Grel[Grel != 2] <- 0
  }
  
  if (graph == TRUE) {
    if(requireNamespace("Rgraphviz", quietly = TRUE, warn.conflicts = FALSE) == TRUE) {
      g <- as( Grel, "graphNEL" )
      plot(g, main = paste("Completed partially directed graph" ) )
    } else {
      warning('In order to plot the generated network, package Rgraphviz is required.')
    }
  }
  list( parents = parents, children = children, spouses = spouses, relatives = relatives, markov.blanket = blanket )
}