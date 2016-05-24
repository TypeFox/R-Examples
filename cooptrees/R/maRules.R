#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences                                   #
#-----------------------------------------------------------------------------#

# maRules ---------------------------------------------------------------------
#' Allocation rules for minimum cost arborescence problems
#' 
#' Given a graph with at least one minimum cost arborescence, the 
#' \code{maRules} function divides the cost of the arborescence among the
#' agents according to "Bird" and "ERO" rules.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param rule denotes the chosen allocation rule: "Bird" or "ERO".
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default its value is
#' \code{TRUE}.
#'
#' @return \code{maRules} returns a matrix with the agents and their costs.
#' It also prints the result in console.
#' 
#' @references B. Dutta and D. Mishra, "Minimum cost arborescences", Games and
#' Economic Behavior, vol. 74, pp. 120-143, Jan. 2012.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,7, 1,3,6, 1,4,4, 2,3,8, 2,4,6, 3,2,6,
#'                  3,4,5, 4,2,5, 4,3,7), ncol = 3, byrow = TRUE)
#' # Allocation rules
#' maRules(nodes, arcs, rule = "Bird")
#' maRules(nodes, arcs, rule = "ERO")

maRules <- function(nodes, arcs, rule, show.data = TRUE) {
  
  # Compute costs according to selected rule
  if (rule == "Bird") {
    coopRules <- maBird(nodes, arcs)
  } else if (rule == "ERO") {
    coopRules <- maERO(nodes, arcs)
  } else {
    warning("Unknown rule")
  }
  
  if (show.data) {
    
    # Build output
    df <- as.data.frame(coopRules)
    df[, 2] <- round(df[, 2], 2)  # round costs
    colnames(df) <- c("Agent ", "  Cost ")
    showdf <- capture.output(print(df, row.names = FALSE))
    pastedf <- paste(showdf, "\n", sep="")
    
    # Show console output
    cat("\n")
    cat(" Minimum cost arborescence \n")
    cat(" Rule:", rule, "\n")
    cat(" ----------------\n")
    cat(" ", pastedf)
    cat(" ----------------\n")
    cat("      Total =", sum(coopRules[, 2]), "\n")
    cat("\n")
    
  }
  
  invisible(coopRules)
  
}
#-----------------------------------------------------------------------------#