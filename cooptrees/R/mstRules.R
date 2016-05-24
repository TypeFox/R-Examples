#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum spanning trees                                       #
#-----------------------------------------------------------------------------#

# mstRules --------------------------------------------------------------------
#' Allocation rules for minimum cost spanning tree problem
#' 
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstRules} function divides the cost of the tree among the agents
#' according to the most known rules: "Bird", "Dutta-Kar", "Kar", "ERO".
#'
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param rule denotes the chosen allocation rule: "Bird", "Dutta-Kar", "Kar"
#' or "ERO".
#' @param algorithm denotes the algorithm used with the ERO rule: "Kruskal".
#' @param show.data logical value indicating if the function displays the 
#' console output\code{TRUE} or no \code{FALSE}. The deafult is \code{TRUE}.
#'
#' @return \code{mstRules} returns a matrix with the agents and the cost
#' that each one of them has to pay. It also prints the result in console.
#' 
#' @references C. G. Bird, "On Cost Allocation for a Spanning Tree: A Game
#' Theoretic Approach", Networks, no. 6, pp. 335-350, 1976.
#' 
#' B. Dutta and A. Kar, "Cost monotonicity, consistency and minimum
#' cost spanning tree games", Games and Economic Behavior, vol. 48,
#' pp. 223-248, Aug. 2004.
#' 
#' A. Kar, "Axiomatization of the Shapley Value on Minimum Cost
#' Spanning Tree Games", Games and Economic Behavior, vol. 38, pp. 265-277,
#' Feb. 2002.
#' 
#' G. Bergantiños and J. J. Vidal-Puga, "A fair rule in minimum
#' cost spanning tree problems", Journal of Economic Theory, vol. 137,
#' pp. 326-352, Nov. 2007.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Allocation Rules
#' mstRules(nodes, arcs, rule = "Bird")
#' mstRules(nodes, arcs, rule = "Dutta-Kar")
#' mstRules(nodes, arcs, rule = "Kar")
#' mstRules(nodes, arcs, rule = "ERO", algorithm = "Kruskal")

mstRules <- function(nodes, arcs, rule, algorithm = "Kruskal",
                     show.data = TRUE) {
  
  # Compute rule according to the option selected
  if (rule == "Bird") {
    coopRules <- mstBird(nodes, arcs)
  } else if (rule == "Dutta-Kar") {
    coopRules <- mstDuttaKar(nodes, arcs)
  } else if (rule == "Kar") {
    coopRules <- mstKar(nodes, arcs)
  } else if (rule == "ERO") {
    if (algorithm == "Kruskal") {
      coopRules <- mstEROKruskal(nodes, arcs)
    } else {
      warning("Unknown algorithm")
    }
  } else {
    warning("Unknown rule")
  }
  
  if (show.data) {
  
    # Output
    df <- as.data.frame(coopRules)
    df[, 2] <- round(df[, 2], 2)  # round costs
    colnames(df) <- c("Agent ", "  Cost ")
    showdf <- capture.output(print(df, row.names = FALSE))
    pastedf <- paste(showdf, "\n", sep="")
    
    # Console output
    cat("\n")
    cat(" Minimum spanning tree \n")
    cat(" Rule:", rule, "\n")
    cat(" ----------------\n")
    cat(" ", pastedf)
    cat(" ----------------\n")
    cat("      Total =", sum(coopRules[, 2]), "\n")
    cat("\n")
  
  }
  
  return(coopRules)
  
}
#-----------------------------------------------------------------------------#