#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost arborescences problems                          #
#-----------------------------------------------------------------------------#

# maCooperative --------------------------------------------------------------
#' Cooperation in minimum cost arborescence problems
#' 
#' Given a graph with at least one minimum cost arborescence, the
#' \code{maCooperative} function computes the cooperative and "Bird" and
#'"ERO" rules.
#' 
#' @param nodes vector containing the nodes of the graph, identified by a 
#' number that goes from \eqn{1} to the order of the graph.
#' @param arcs matrix with the list of arcs of the graph. Each row represents
#' one arc. The first two columns contain the two endpoints of each arc and the
#' third column contains their weights.
#' @param show.data logical value indicating if the function displays the 
#' console output (\code{TRUE}) or not (\code{FALSE}). By default its value is
#' \code{TRUE}.
#'
#' @return \code{maCooperative} returns and prints a list with the cooperative
#' games and allocation rules of a minimum cost arborescence problem.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,2, 1,3,3, 1,4,4, 2,3,3, 2,4,4, 3,2,3, 
#'                  3,4,1, 4,2,1, 4,3,2), ncol = 3, byrow = TRUE)
#' # Cooperation in minimum cost arborescence problems
#' maCooperative(nodes, arcs)

maCooperative <- function(nodes, arcs, show.data = TRUE) {
  
  # Number of players
  n <- length(nodes) - 1
  
  # Cooperative games
  pesGame <- maGames(nodes, arcs, show.data = FALSE)
  # Irreducible form
  maIrr <- maIrreducible(nodes, arcs)
  pesIrrGame <- maGames(nodes, maIrr, show.data = FALSE)
  
  # Rules
  # Bird
  bRule <- maRules(nodes, arcs, rule = "Bird", show.data = FALSE)
  bRule <- bRule[order(bRule[, 1]), ]
  # ERO
  ekRule <- maRules(nodes, arcs, rule = "ERO", show.data = FALSE)
  ekRule <- ekRule[order(ekRule[, 1]), ]
  
  # If there is 3 players make a graph with games and rules
  if (n == 3) {
    # Cost game
    v <- -pesGame$values
    # Imputation set
    imputations <- getImputationSet(3, v)
    plotImputationSet(imputations)
    # Core set
    core <- getCoreSet(3, v)
    plotCoreSet(core)
    
    # Rules
    # Bird
    bRule2D <- point3Dto2D(-bRule[, 2], 30, 30)
    points(bRule2D[1], bRule2D[2], pch = 15, col = "yellow", cex = 1.5)
    # ERO  
    ekRule2D <- point3Dto2D(-ekRule[, 2], 30, 30)
    points(ekRule2D[1], ekRule2D[2], pch = 18, col = "red", cex = 1.7)
    
    # Legend
    legend("topright", c("Bird", "ERO"), cex = 0.8, bty = "n",
           col = c("yellow", "red"), pch = c(15, 18), pt.cex = 1.5,
           y.intersp = 1.25, inset = c(0.05, 0.025))
  }
  
  if (show.data) {
  
    # Cooperative games matrix
    gamesMat <- matrix(c(pesGame$values, pesIrrGame$values),
                       nrow = 2, byrow = TRUE)
    # Output
    df1 <- as.data.frame(gamesMat)
    colnames(df1) <- pesGame$coalitions
    rownames(df1) <- c("Pes.", "Irr.")
    
    # Allocation rules matrix
    rulesMat <- matrix(c(bRule[, 2], ekRule[, 2]),
                       nrow = n)
    # Output
    df <- as.data.frame(rulesMat)
    df[, 1:2] <- round(df[, 1:2], 2)  # round costs
    colnames(df) <- c("Bird", "ERO")
    rownames(df) <- nodes[-1]
    
    # Console output
    getMinimumArborescence(nodes, arcs, show.graph = FALSE)
    cat(" Cooperative games: \n")
    cat(" ---------------------------\n")
    print(df1)
    cat(" ---------------------------\n")
    cat("\n")
    cat(" Allocation rules: \n")
    cat(" ---------------------------\n")
    print(df)
    cat(" ---------------------------\n")
    cat("   Total cost =", sum(rulesMat[, 2]), "\n")
    cat("\n")

  }
  
  invisible(list("pessimistic" = pesGame$values,
                 "Bird" = bRule[, 2], "ERO" = ekRule[, 2]))
  
}
#-----------------------------------------------------------------------------#