#-----------------------------------------------------------------------------#
# optrees Package                                                             #
# Minimun Arborescence Problems                                               #
#-----------------------------------------------------------------------------#

# getMinimumArborescence ------------------------------------------------------
#' Computes a minimum cost arborescence
#'
#' Given a connected weighted directed graph, \code{getMinimumArborescence} 
#' computes a minimum cost arborescence. This function provides a method 
#' to find the minimum cost arborescence with Edmonds' algorithm.
#' 
#' @details Given a connected weighted directed graph, a minimum cost 
#' arborescence is an arborescence such that the sum of the weight of its arcs
#' is minimum. In some cases, it is possible to find several minimum cost 
#' arborescences, but the proposed algorithm only finds one of them.
#' 
#' Edmonds' algorithm was developed by the mathematician and computer
#' scientist Jack R. Edmonds in 1967. Although, it was previously proposed in
#' 1965 by Yoeng-jin Chu and Tseng-hong Liu. This algorithm decreases the
#' weights of the arcs in a graph and compacts cycles of zero weight until it
#' can find an arborescence. This arborescence has to be a minimum cost
#' arborescence of the graph.
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param source.node number pointing to the source node of the graph. It's node
#' \eqn{1} by default.
#' @param algorithm denotes the algorithm used to find a minimum cost
#' arborescence: "Edmonds".
#' @param check.graph logical value indicating if it is necesary to check the
#' graph. Is \code{FALSE} by default.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). The default is
#' \code{TRUE}.
#' @param show.graph logical value indicating if the function displays a 
#' graphical representation of the graph and its minimum arborescence
#' (\code{TRUE}) or not (\code{FALSE}). The default is \code{TRUE}.
#' @param stages.data logical value indicating if the function returns data of
#' each stage. The default is \code{FALSE}.
#'
#' @return \code{getMinimumArborescence} returns a list with:
#' tree.nodes vector containing the nodes of the minimum cost arborescence.
#' tree.arcs matrix containing the list of arcs of the minimum cost 
#' arborescence.
#' weight value with the sum of weights of the arcs.
#' stages number of stages required.
#' time time needed to find the minimum cost arborescence.
#' This function also represents the graph and the minimum arborescence and
#' prints to the console the results with additional information (number of
#' stages, computational time, etc.).
#' 
#' @references Chu, Y. J., and Liu, T. H., "On the Shortest Arborescence of a
#' Directed Graph", Science Sinica, vol. 14, 1965, pp. 1396-1400.
#' 
#' Edmonds, J., "Optimum Branchings", Journal of Research of the National
#' Bureau of Standards, vol. 71B, No. 4, October-December 1967, pp. 233-240.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,2, 1,3,3, 1,4,4, 2,3,3, 2,4,4, 3,2,3, 
#'                  3,4,1, 4,2,1, 4,3,2),byrow = TRUE, ncol = 3)
#' # Minimum cost arborescence
#' getMinimumArborescence(nodes, arcs)
#' 
#' @export

getMinimumArborescence <- function(nodes, arcs, source.node = 1,
                                   algorithm = "Edmonds", stages.data = FALSE,
                                   show.data = TRUE, show.graph = TRUE,
                                   check.graph = FALSE) {
  
  # Check validity of the graph
  if (check.graph) {
    validGraph <- checkGraph(nodes, arcs, source.node, directed = TRUE)
    if (validGraph == FALSE) {
      # If graph is not valid stop algorithm before start
      validGraph <- "Not valid graph"
      return(validGraph)
    }
  }    
  
  # Remove loops
  arcs <- removeLoops(arcs)
  # Remove multi-arcs
  arcs <- removeMultiArcs(arcs, directed = TRUE)
  
  # Always check time
  tini <- proc.time()  # start time
  
  algorithm <- "Edmonds"
  msa <- msArborEdmonds(nodes, arcs, source.node, stages.data)
  
  running.time <- proc.time() - tini  # elapsed time
  
  # Console output
  if (show.data) {
    # Build data frame to show
    dfArbor <- msa$tree.arcs
    weight <- sum(dfArbor[, 3])
    df <- as.data.frame(dfArbor)
    colnames(df) <- c("   head ", "   tail ", "  weight ")
    showdf <- capture.output(print(df, row.names=FALSE))
    pastedf <- paste(showdf, "\n", sep="")
    
    # Output
    cat("\n")
    cat(" Minimum cost spanning arborescence \n")
    cat(" Algorithm:", algorithm, "\n")
    cat(" Stages: ", msa$stages, "| Time: ", running.time[3], "\n")
    cat(" ------------------------------\n")
    cat(" ", pastedf)
    cat(" ------------------------------\n")
    cat("                   Total =", weight, "\n")
    cat("\n")
    
  }
  
  # Graph
  if (show.graph) {
    repGraph(nodes, arcs, tree = msa$tree.arcs, directed = TRUE,
             plot.title = "Minimum Cost Spanning Arborescence")
  }
  
  # Return
  if (stages.data == TRUE) {
    output <- list("tree.nodes" = msa$tree.nodes, "tree.arcs" = msa$tree.arcs,
                   "weight" = sum(msa$tree.arcs[, 3]), "stages" = msa$stages,
                   "stage" = msa$stages.data, "time" = running.time)
  } else {
    output <- list("tree.nodes" = msa$tree.nodes, "tree.arcs" = msa$tree.arcs,
                   "weight" = sum(msa$tree.arcs[, 3]),
                   "stages" = msa$stages, "time" = running.time)
  }
  
  invisible(output)
  
}
#-----------------------------------------------------------------------------#