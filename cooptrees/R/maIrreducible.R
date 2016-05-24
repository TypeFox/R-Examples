#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences problems                          #
#-----------------------------------------------------------------------------#

# maIrreducible ---------------------------------------------------------------
#' Irreducible form for a minimum cost arborescence problem
#' 
#' Given a graph with at least one minimum cost arborescence the
#' \code{maIrreducible} function obtains the irreducible form.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#'
#' @return \code{maIrreducible} returns a matrix with the list of arcs of the
#' irreducible form.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
#' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # Irreducible form          
#' maIrreducible(nodes, arcs)

maIrreducible <- function(nodes, arcs) {
  
  # Get minimum cost arborescence and every step data
  mca <- getMinimumArborescence(nodes, arcs, stages.data = TRUE,
                                show.data = FALSE, show.graph = FALSE)
  
  # Number of steps to find it
  numStages <- mca$stages
  i <- numStages; i
  
  # Partir de la etapa inicial
  prevArcs <- mca$stage[[i]]$arcs
  # Extract mcArcs
  mcArcs <- mca$stage[[i]]$mcArcs
  # Cada arco se queda con el coste mínimo que llega a su node de destino
  for (k in 1:nrow(mcArcs)) {
    prevArcs[which(prevArcs[, 2] == mcArcs[k, 2]), 3] <- mcArcs[k, 3]
  }
  # Guardamos prevArcs por si hay solo una etapa
  stageArcs <- prevArcs
  
  i <- i - 1
  # Iteramos etapa por etapa
  while (i > 0) {
    
    # Arcos de la nueva etapa
    stageArcs <- mca$stage[[i]]$arcs
    # Supernodo
    superNode <- mca$stage[[i]]$superNode
    # Correspondencias  
    matchesNodes <- mca$stage[[i]]$matches  # correspondences between nodes
    # Empezar con coste 0
    stageArcs[, 3] <- 0
    
    # Revisamos arcos previos
    for (k in 1:nrow(prevArcs)) {
      # Si el nodo de partida o el nodo de llegada son un supernodo
      if (prevArcs[k, 1] %in% superNode) {
        # Qué arcos de los nuevos se corresponden
        # Nodos de partida
        leavingNodes <- matchesNodes[which(matchesNodes[, 2] == superNode), 1]
        # Correspondencia nodos de llegada
        reachingNodes <- matchesNodes[matchesNodes[, 2] == prevArcs[k, 2], 1]
        # Arcos que tengan los nodos de partida y el de llegada
        matchesArcs <- which(stageArcs[, 2] %in% reachingNodes 
                             & stageArcs[, 1] %in% leavingNodes)
        # A esos arcos se les cambia el peso por el del arco previo
        stageArcs[matchesArcs, 3] <- prevArcs[k, 3]    
      } else if (prevArcs[k, 2] %in% superNode) {
        # Qué arcos de los nuevos se corresponden
        # Nodos de llegada
        reachingNodes <- matchesNodes[which(matchesNodes[, 2] == superNode), 1]
        # Correspondencia nodos de partida
        leavingNodes <- matchesNodes[matchesNodes[, 2] == prevArcs[k, 1], 1]
        # Arcos que tengan el nodo de partida y los de llegada
        matchesArcs <- which(stageArcs[, 1] %in% leavingNodes 
                             & stageArcs[, 2] %in% reachingNodes)
        # A esos arcos se les cambia el peso por el del arco previo
        stageArcs[matchesArcs, 3] <- prevArcs[k, 3]
      } else {
        # Nodos sin relación con el supernodo
        # Correspondencia nodos de partida
        leavingNodes <- matchesNodes[matchesNodes[, 2] == prevArcs[k, 1], 1]
        # Correspondencia nodos de llegada
        reachingNodes <- matchesNodes[matchesNodes[, 2] == prevArcs[k, 2], 1]
        stageArcs[stageArcs[, 1] == leavingNodes &
                    stageArcs[, 2] == reachingNodes, 3] <- prevArcs[k, 3]
      }
    }
    
    # Extract mcArcs
    mcArcs <- mca$stage[[i]]$mcArcs
    # Nodos de destino
    destNodes <- unique(mcArcs[, 2])
    # A cada arco se le suma el coste mínimo que llega a su nodo de destino
    for (k in 1:length(destNodes)) {
      prevCost <- stageArcs[which(stageArcs[, 2] == destNodes[k]), 3]
      newCost <- mcArcs[which(mcArcs[, 2] == destNodes[k])[1], 3]
      stageArcs[which(stageArcs[, 2] == destNodes[k]), 3] <- prevCost + newCost
    }
    
    # Los arcos de la etapa se convierten en arcos previos
    prevArcs <- stageArcs
    # Y repetimos
    i <- i - 1
  }
  
  # Finish with list of arcs of the irreducible form
  irreducibleForm <- stageArcs
  
  # Return irreducible form
  return(irreducibleForm)
  
}
#-----------------------------------------------------------------------------#