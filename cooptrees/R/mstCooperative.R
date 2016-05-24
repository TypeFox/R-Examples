#-----------------------------------------------------------------------------#
# Cooptrees package                                                           #
# Cooperation in minimum cost spanning trees problems                         #
#-----------------------------------------------------------------------------#

# mstCooperative --------------------------------------------------------------
#' Cooperation in minimum cost spanning tree problems
#' 
#' Given a graph with at least one minimum cost spanning tree, the
#' \code{mstCooperative} function computes the pessimistic and optimistic
#' games; and the most known allocation rules: "Bird", "Dutta-Kar","Kar" and
#' "ERO".
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
#' @return \code{mstCooperative} returns and print a list with the cooperative
#' games and the allocation rules of a minimum cost spanning tree problem.
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
#' G. Bergantiños and J. J. Vidal-Puga, "The optimistic TU game in
#' minimum cost spanning tree problems", International Journal of Game Theory,
#' vol. 36, pp. 223-239, Feb. 2007.
#' 
#' V. Feltkamp, S. H. Tijs, S. Muto, "On the irreducible core and the
#' equal remaining obligation rule of minimum cost extension problems",
#' Mimeo, Tilburg University, 1994.
#' 
#' @examples
#' # Graph
#' nodes <- 1:4
#' arcs <- matrix(c(1,2,6, 1,3,10, 1,4,6, 2,3,4, 2,4,6, 3,4,4), 
#'                byrow = TRUE, ncol = 3)
#' # Cooperation in minimum cost spanning tree problems
#' mstCooperative(nodes, arcs)

mstCooperative <- function(nodes, arcs, show.data = TRUE) {
  
  # Number of players
  n <- length(nodes) - 1
  
  # Cooperative games
  pesGame <- mstGames(nodes, arcs, game = "pessimistic",
                      show.data = FALSE)
  optGame <- mstGames(nodes, arcs, game = "optimistic",
                      show.data = FALSE)
  # Irreducible form
  mstIrr <- mstIrreducible(nodes, arcs)
  pesIrrGame <- mstGames(nodes, mstIrr, game = "pessimistic",
                         show.data = FALSE)
  
  # Rules
  # Bird
  bRule <- mstRules(nodes, arcs, rule = "Bird", show.data = FALSE)
  bRule <- bRule[order(bRule[, 1]), ]
  # Dutta-Kar
  dkRule <- mstRules(nodes, arcs, rule = "Dutta-Kar", show.data = FALSE)
  dkRule <- dkRule[order(dkRule[, 1]), ]
  # Kar
  kRule <- mstRules(nodes, arcs, rule = "Kar", show.data = FALSE)
  kRule <- kRule[order(kRule[, 1]), ]
  # ERO
  ekRule <- mstRules(nodes, arcs, rule = "ERO", show.data = FALSE)
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
    # Dutta-Kar
    dkRule2D <- point3Dto2D(-dkRule[, 2], 30, 30)
    points(dkRule2D[1], dkRule2D[2], pch = 16, col = "green", cex = 1.5)
    # Kar
    kRule2D <- point3Dto2D(-kRule[, 2], 30, 30)
    points(kRule2D[1], kRule2D[2], pch = 17, col = "blue", cex = 1.5)
    # ERO
    ekRule2D <- point3Dto2D(-ekRule[, 2], 30, 30)
    points(ekRule2D[1], ekRule2D[2], pch = 18, col = "red", cex = 1.7)
    
    # Legend
    legend("topright", c("Bird", "Dutta-Kar", "Kar", "ERO"), bty = "n",
           col = c("yellow", "green", "blue", "red"), pch = c(15, 16, 17, 18), 
           pt.cex = 1.5, y.intersp = 1.25, inset = c(0.05, 0.025), cex = 0.8)
  }
  
  if (show.data) {
  
    # Cooperative games matrix
    gamesMat <- matrix(c(pesGame$values, optGame$values, pesIrrGame$values),
                       nrow = 3, byrow = TRUE)
    # Output
    df1 <- as.data.frame(gamesMat)
    colnames(df1) <- pesGame$coalitions
    rownames(df1) <- c("Pes.", "Opt.", "Irr.")
    
    # Allocation rules matrix
    rulesMat <- matrix(c(bRule[, 2], dkRule[, 2], kRule[, 2], ekRule[, 2]),
                       nrow = n)
    # Output
    df <- as.data.frame(rulesMat)
    df[, 1:4] <- round(df[, 1:4], 2)  # round costs
    colnames(df) <- c("Bird", "D-K", "Kar", "ERO")
    rownames(df) <- nodes[-1]
    
    # Console output
    getMinimumSpanningTree(nodes, arcs, algorithm = "Prim", show.graph = FALSE)
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

  invisible(list("pessimistic" = pesGame$values, "optimistic" = optGame$values,
              "Bird" = bRule[, 2], "Dutta-Kar" = dkRule[, 2],
              "Kar" = kRule[, 2], "ERO" = ekRule[, 2]))
  
}
#-----------------------------------------------------------------------------#