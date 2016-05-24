#-----------------------------------------------------------------------------#
# Optrees Package                                                             #
# Auxiliar functions                                                          #
#-----------------------------------------------------------------------------#

# repGraph --------------------------------------------------------------------
#' Visual representation of a graph
#' 
#' The \code{repGraph} function uses \code{igraph} package to represent a 
#' graph. 
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param tree matrix with the list of arcs of a tree, if there is one. Is 
#' \code{NULL} by default.
#' @param directed logical value indicating whether the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' @param plot.title string with main title of the graph. Is \code{NULL} by
#' default.
#' @param fix.seed number to set a seed for the representation.
#' 
#' @return \code{repGraph} returns a plot with the given graph.
#' 
#' @examples
#' # Graph
#' nodes <- c(1:4)
#' arcs <- matrix(c(1,2,2, 1,3,15, 2,3,1, 2,4,9, 3,4,1), 
#'                byrow = TRUE, ncol = 3)
#' # Plot graph
#' repGraph(nodes, arcs)
#' 
#' @export

repGraph <- function(nodes, arcs, tree = NULL, directed = FALSE,
                     plot.title = NULL, fix.seed = NULL) {
  
  #require(igraph)  # need igraph library
  ref.size <- 30 - length(nodes)

  g <- graph(t(arcs[, 1:2]), directed = directed)  # create igraph object
  E(g)$weight <- arcs[, 3]  # add weights

  # Format
  V(g)$size <- ref.size # size of the nodes
  V(g)$color <- "white"  # color of the nodes
  V(g)$label.cex <- 1.2 - length(nodes)*0.03
  if (length(nodes) >= 20) {
    vLabel <- NA
  } else {
    vLabel <- nodes
  }
  # Edges
  eLabelcex <- 1 - nrow(arcs)*0.01
  if (eLabelcex >= 0.5) {
    eLabel <- E(g)$weight
  } else {
    eLabel <- NA
  }
  
  if (!is.null(tree)) {
    
    if (!directed) {
      if (nrow(arcs) == 1) {
        # If there is only one arc
        E(g)$color <- "red"
        E(g)$lty <- 1
      } else {
        # Mark arcs of the tree
        arcList1 <- apply(matrix(arcs[, 1:2], ncol = 2), 1, paste, collapse="-")
        arcList2 <- apply(matrix(tree[, 1:2], ncol = 2), 1, paste, collapse="-")
        arcList3 <- apply(matrix(tree[, 2:1], ncol = 2), 1, paste, collapse="-")
        E(g)$color <- ifelse(arcList1 %in% c(arcList2, arcList3), "red", "gray")
        E(g)$lty <- ifelse(arcList1 %in% c(arcList2, arcList3), 1, 1)
      }
    
    } else {
      
      if (nrow(arcs) == 1) {
        # If there is only one arc
        E(g)$color <- "red"
        E(g)$lty <- 1
      } else {
        # Mark arcs of the msTree
        arcList1 <- apply(matrix(arcs[, 1:2], ncol = 2), 1, paste, collapse="-")
        if (nrow(tree) == 1) {
          arcList2 <- paste(matrix(tree[, 1:2], ncol = 2), collapse = "-")
        } else {
          arcList2 <- apply(matrix(tree[, 1:2], ncol = 2), 1, paste, collapse="-")
        }
        E(g)$color <- ifelse(arcList1 %in% arcList2, "red", "gray")
        E(g)$lty <- ifelse(arcList1 %in% arcList2, 1, 1)
      }
      
    }
    
  }

  # Seed
  if (!is.null(fix.seed)) {
    set.seed(fix.seed)
  }
  
  # Plot complete graph
  plot.igraph(g, frame = FALSE, main = plot.title, vertex.label = vLabel,
              vertex.label.color = "black", vertex.label.family = "sans",
              edge.curved = directed, edge.label = eLabel, edge.width = 2,
              edge.label.color="black", edge.label.cex=eLabelcex,
              edge.label.family = "sans", edge.arrow.size = 0.3,
              edge.label.font = 1)
  
}
#-----------------------------------------------------------------------------#