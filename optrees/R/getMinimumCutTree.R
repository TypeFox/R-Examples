#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimum Cut Tree Problems                                                   #
#-----------------------------------------------------------------------------#

#' getMinimumCutTree ----------------------------------------------------------
#' Computes a minimum cut tree
#'
#' Given a connected weighted undirected graph, \code{getMinimumCutTree}
#' computes a minimum cut tree, also called Gomory-Hu tree. This function
#' uses the Gusfield's algorithm to find it.
#' 
#' @details The minimum cut tree or Gomory-Hu tree was introduced by R. E. 
#' Gomory and T. C. Hu in 1961. Given a connected weighted undirected graph, 
#' the Gomory-Hu tree is a weighted tree that contains the minimum s-t cuts 
#' for all s-t pairs of nodes in the graph. Gomory and Hu developed an 
#' algorithm to find this tree, but it involves maximum flow searchs and nodes
#' contractions.
#' 
#' In 1990, Dan Gusfield proposed a new algorithm that can be used to find the
#' Gomory-Hu tree without any nodes contraction and simplifies the 
#' implementation.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param algorithm denotes the algorithm to use for find a minimum cut tree or
#' Gomory-Hu tree: "Gusfield".
#' @param check.graph logical value indicating if it is necesary to check the
#' graph. Is \code{FALSE} by default.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). The default is
#' \code{TRUE}.
#' @param show.graph logical value indicating if the function displays a 
#' graphical representation of the graph and its minimum cut tree
#' (\code{TRUE}) or not (\code{FALSE}). The default is \code{TRUE}.
#' 
#' @return \code{getMinimumCutTree} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cut tree.
#' tree.arcs matrix containing the list of arcs of the minimum cut tree.
#' weight value with the sum of weights of the arcs.
#' stages number of stages required.
#' time time needed to find the minimum cost spanning tree.
#' This function also represents the graph and the minimum cut tree and
#' prints in console the results whit additional information (number of stages,
#' computational time, etc.).
#' 
#' @references R. E. Gomory, T. C. Hu. Multi-terminal network flows. Journal 
#' of the Society for Industrial and Applied Mathematics, vol. 9, 1961.
#' 
#' Dan Gusfield (1990). "Very Simple Methods for All Pairs Network Flow 
#' Analysis". SIAM J. Comput. 19 (1): 143-155.
#' 
#' @examples
#' # Graph
#' nodes <- 1:6
#' arcs <- matrix(c(1,2,1, 1,3,7, 2,3,1, 2,4,3, 2,5,2, 3,5,4, 4,5,1, 4,6,6, 
#'                 5,6,2), byrow = TRUE, ncol = 3)
#' # Minimum cut tree
#' getMinimumCutTree(nodes, arcs)
#' 
#' @export

getMinimumCutTree <- function(nodes, arcs, algorithm = "Gusfield",
                              show.data = TRUE, show.graph = TRUE,
                              check.graph = FALSE) {
                     
  # Check validity of the graph
  if (check.graph) {
    validGraph <- checkGraph(nodes, arcs, source.node = 1, directed = FALSE)
    if (validGraph == FALSE) {
      # If graph is not valid stop algorithm before start
      validGraph <- "Invalid graph"
      return(validGraph)
    }
  }
  
  # Remove loops
  arcs <- removeLoops(arcs)
  # Remove multi-arcs
  arcs <- removeMultiArcs(arcs, directed = FALSE)
  
  # Always check time
  tini <- proc.time()  # start time
  
  if (algorithm == "Gusfield") {
    ght <- ghTreeGusfield(nodes, arcs)
  } else {
    stop("Unknown algorithm")
  }
  
  running.time <- proc.time() - tini  # elapsed time
  
  # Console output
  if (show.data) {
    # Build data frame to show
    dfTree <- ght$tree.arcs
    weight <- sum(dfTree[, 3])
    df <- as.data.frame(dfTree)
    colnames(df) <- c("   ept1 ", "   ept2 ", "  weight ")
    showdf <- capture.output(print(df, row.names=FALSE))
    pastedf <- paste(showdf, "\n", sep="")
    
    # Output
    cat("\n")
    cat(" Minimum cut tree (Gomory-Hu) \n")
    cat(" Algorithm:", algorithm, "\n")
    cat(" Steps: ", ght$stages, "| Time: ", running.time[3], "\n")
    cat(" ------------------------------\n")
    cat(" ", pastedf)
    cat(" ------------------------------\n")
    cat("                    Total =", weight, "\n")
    cat("\n")
    
  }
  
  # Graph
  if (show.graph) {
    
    # Show two graphs
    par(mfrow = c(1, 2))
    
    # Original graph
    repGraph(nodes, arcs, plot.title = "Original Graph", fix.seed = 1)
    # Minimum cut tree
    repGraph(nodes, ght$tree.arcs, ght$tree.arcs,
             plot.title = "Minimum Cut Tree", fix.seed = 1)
    
    # Restore default R graphic parameters
    par(mfrow = c(1, 1))
    
  }
  
  # Return
  output <- list("tree.nodes" = ght$tree.nodes, "tree.arcs" = ght$tree.arcs,
                 "weight" = sum(ght$tree.arcs[, 3]),
                 "stages" = ght$stages, "time" = running.time)
  invisible(output)
  
}
#-----------------------------------------------------------------------------#