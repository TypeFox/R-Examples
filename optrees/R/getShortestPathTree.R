#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Shortest Path Tree Problems                                                 #
#-----------------------------------------------------------------------------#

# getShortestPathTree ---------------------------------------------------------
#' Computes a shortest path tree
#'
#' Given a connected weighted graph, directed or not,
#' \code{getShortestPathTree} computes the shortest path tree from a given 
#' source node to the rest of the nodes the graph, forming a shortest path
#' tree. This function provides methods to find it with two known algorithms:
#' "Dijkstra" and "Bellman-Ford".
#' 
#' @details Given a connected weighted graph, directed or not, a shortest
#' path tree rooted at a source node is a spanning tree such thtat the path
#' distance from the source to any other node is the shortest path distance
#' between them. Differents algorithms were proposed to find a shortest path
#' tree in a graph.
#' 
#' One of these algorithms is Dijkstra's algorithm. Developed by the computer
#' scientist Edsger Dijkstra in 1956 and published in 1959, it is an algorithm
#' that can compute a shortest path tree from a given source node to the others
#' nodes in a connected, directed or not, graph with non-negative 
#' weights.
#' 
#' The Bellman-Ford algorithm gets its name for two of the developers,
#' Richard Bellman y Lester Ford Jr., and it was published by them in 1958 and
#' 1956 respectively. The same algorithm also was published independently in
#' 1957 by Edward F. Moore. This algorithm can compute the shortest path from a
#' source node to the rest of nodes in a connected, directed or not,
#' graph with weights that can be negatives. If the graph is connected and
#' there isn't negative cycles, the algorithm always finds a shortest path
#' tree.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param algorithm denotes the algorithm used to find a shortest path tree:
#' "Dijkstra" or "Bellman-Ford".
#' @param check.graph logical value indicating if it is necesary to check the
#' graph. Is \code{FALSE} by default.
#' @param source.node number indicating the source node of the graph. It's node
#' \eqn{1} by default.
#' @param directed logical value indicating wheter the graph is directed 
#' (\code{TRUE}) or not (\code{FALSE}).
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). The default is
#' \code{TRUE}.
#' @param show.graph logical value indicating if the function displays a 
#' graphical representation of the graph and its shortest path tree
#' (\code{TRUE}) or not (\code{FALSE}). The default is \code{TRUE}.
#' @param show.distances logical value indicating if the function displays in
#' the console output the distances from source to all the other nodes. The
#' default is \code{TRUE}.
#'
#' @return \code{getShortestPathTree} returns a list with:
#' tree.nodes vector containing the nodes of the shortest path tree.
#' tree.arcs matrix containing the list of arcs of the shortest path tree.
#' weight value with the sum of weights of the arcs.
#' distances vector with distances from source to other nodes
#' stages number of stages required.
#' time time needed to find the minimum cost spanning tree.
#' This function also represents the graph with the shortest path tree 
#' and prints to the console the results with additional information (number of
#' stages, computational time, etc.).
#' 
#' @references Dijkstra, E. W. (1959). "A note on two problems in connexion 
#' with graphs". Numerische Mathematik 1, 269-271.
#' 
#' Bellman, Richard (1958). "On a routing problem". Quarterly of 
#' Applied Mathematics 16, 87-90.
#' 
#' Ford Jr., Lester R. (1956). Network Flow Theory. Paper P-923. Santa Monica, 
#' California: RAND Corporation.
#' 
#' Moore, Edward F. (1959). "The shortest path through a maze". Proc. Internat.
#' Sympos. Switching Theory 1957, Part II. Cambridge, Mass.: Harvard Univ.
#' Press. pp. 285-292.
#' 
#' @examples
#' # Graph
#' nodes <- 1:5
#' arcs <- matrix(c(1,2,2, 1,3,2, 1,4,3, 2,5,5, 3,2,4, 3,5,3, 4,3,1, 4,5,0),
#'                ncol = 3, byrow = TRUE)
#' # Shortest path tree
#' getShortestPathTree(nodes, arcs, algorithm = "Dijkstra", directed=FALSE)
#' getShortestPathTree(nodes, arcs, algorithm = "Bellman-Ford", directed=FALSE)
#' 
#' # Graph with negative weights
#' nodes <- 1:5
#' arcs <- matrix(c(1,2,6, 1,3,7, 2,3,8, 2,4,5, 2,5,-4, 3,4,-3, 3,5,9, 4,2,-2,
#'                  5,1,2, 5,4,7), ncol = 3, byrow = TRUE)
#' # Shortest path tree
#' getShortestPathTree(nodes, arcs, algorithm = "Bellman-Ford", directed=TRUE)
#' 
#' @export

getShortestPathTree <- function(nodes, arcs, algorithm, check.graph = FALSE,
                                source.node = 1, directed = TRUE,
                                show.data = TRUE, show.graph = TRUE,
                                show.distances = TRUE) {
  
  # Check validity of the graph
  if (check.graph) {
    validGraph <- checkGraph(nodes, arcs, source.node, directed)
    if (validGraph == FALSE) {
      # If graph is not valid stop algorithm before start
      validGraph <- "Not a valid graph"
      return(validGraph)
    }
  }
  
  # Remove loops
  arcs <- removeLoops(arcs)
  # Remove multi-arcs
  arcs <- removeMultiArcs(arcs, directed)
  
  # Always check time
  tini <- proc.time()  # start time
  
  if (algorithm == "Dijkstra") {
    spt <- spTreeDijkstra(nodes, arcs, source.node, directed)
  } else if (algorithm == "Bellman-Ford") {
    spt <- spTreeBellmanFord(nodes, arcs, source.node, directed)
  } else {
    stop("Unknown algorithm")
  }
  
  running.time <- proc.time() - tini  # elapsed time  
  
  # Console output
  if (show.data) {
    
    # Build data frame to show
    dfTree <- spt$tree.arcs
    weight <- sum(dfTree[, 3])
    df <- as.data.frame(dfTree)
    if (directed) {
      colnames(df) <- c("   head ", "   tail ", "  weight ")
    } else {
      colnames(df) <- c("   ept1 ", "   ept2 ", "  weight ")
    }
    showdf <- capture.output(print(df, row.names=FALSE))
    pastedf <- paste(showdf, "\n", sep="")
    
    # Output
    cat("\n")
    cat(" Shortest path tree \n")
    cat(" Algorithm:", algorithm, "\n")
    cat(" Stages: ", spt$stages, "| Time: ", running.time[3], "\n")
    cat(" ------------------------------\n")
    cat(" ", pastedf)
    cat(" ------------------------------\n")
    cat("                   Total =", weight, "\n")
    cat("\n")
    
  }
  
  if (show.distances) {
    
    # Build data frame to show
    dfDist <- matrix(c(rep(source.node, length(nodes[-source.node])),
                       nodes[-which(nodes == source.node)],
                       spt$distances[-source.node]), ncol = 3)
    df2 <- as.data.frame(dfDist)
    colnames(df2) <- c(" source ", "   node ", "distance ")
    showdf2 <- capture.output(print(df2, row.names=FALSE))
    pastedf2 <- paste(showdf2, "\n", sep="")
    #cat("\n")
    cat(" Distances from source: \n")
    cat(" ------------------------------\n")
    cat(" ", pastedf2)
    cat(" ------------------------------\n")
    cat("\n")
    
  }
  
  # Graph
  if (show.graph) {
    repGraph(nodes, arcs, tree = spt$tree.arcs, directed = directed,
             plot.title = "Shortest Path Tree")
  }
  
  # Return
  output <- list("tree.nodes" = spt$tree.nodes, "tree.arcs" = spt$tree.arcs,
            "weight" = sum(spt$tree.arcs[, 3]), "distances" = spt$distances,
                 "stages" = spt$stages, "time" = running.time)
  invisible(output)
  
}
#-----------------------------------------------------------------------------#