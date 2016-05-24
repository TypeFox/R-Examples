#-----------------------------------------------------------------------------#
# cooptrees Package                                                           #
# Helper functions                                                            #
#-----------------------------------------------------------------------------#

# shapleyValue ----------------------------------------------------------------
#' Shapley value of a cooperative game
#'
#' Given a cooperative game, the \code{shapleyValue} function computes its 
#' Shapley value.
#'
#' The Shapley value is a solution concept in cooperative game theory proposed
#' by Lloyd Shapley in 1953. It is obtained as the average of the marginal
#' contributions of the players associated with all the posible orders
#' of the players.
#'
#' @param n number of players in the cooperative game.
#' @param S vector with all the possible coalitions. If none has been
#' specified the function generates it automatically.
#' @param v vector with the characteristic function of the cooperative game.
#'
#' @return The \code{shapleyValue} functions returns a matrix with all the
#' marginal contributions of the players (\code{contributions}) and a vector
#' with the Shapley value (\code{value}).
#' 
#' @references Lloyd S. Shapley. "A Value for n-person Games". In Contributions
#' to the Theory of Games, volume II, by H.W. Kuhn and A.W. Tucker, editors.
#' Annals of Mathematical Studies v. 28, pp. 307-317. Princeton University
#' Press, 1953.
#' 
#' @examples
#' # Cooperative game
#' n <- 3  # players
#' v <- c(4, 4, 4, 8, 8, 8, 14)  # characteristic function
#' # Shapley value
#' shapleyValue(n, v = v)

shapleyValue <- function(n, S = NULL, v) {
  
  # Alert when the number of players is big
  continue <- "Y"
  if (n > 9) {
    cat("shapley need some time to compute a", n, "player game \n")
    cat("Continue? (Y = Yes, N = No) \n")
    continue <- scan(what = "character", n = 1)
  }
  
  # Proceed only if it has been confirmed
  if (continue == "Y") {
      
    # List of players in the cooperative game
    players <- 1:n
    
    # Generate coalitions if have not been given
    if (is.null(S)) {
      S <- c()
      for (i in players) {  # compute all possible combinations of i players
        coalitions <- combn(length(players), i)
        for (j in 1:ncol(coalitions)) {
          S <- c(S, paste(coalitions[, j], sep = "", collapse = ","))
        }
      }
    }
    
    # Generate all possible orders of the n players
    orders <- permutations(n = n, r = n, v = players)
    
    # Contributions matrix to store calculations for each player and order
    contributions <- matrix(nrow = nrow(orders), ncol = n)
    
    # Computes values for each order
    for (i in 1:nrow(orders)) {
      
      # Extract row with one order
      actualOrder <- orders[i, ]
      # Vector to store each coalition according to its order
      coalPlayers <- c()
      
      # Do math for each formed coalition according to the order
      for (j in 1:length(actualOrder)) {
        
        # Extract players that forms the coalition
        coalPlayers <- c(coalPlayers, actualOrder[j])
        coal <- paste(coalPlayers[order(coalPlayers)], sep = "", collapse = ",")
        # Select corresponding coalition
        k <- which(S == coal)
        # Subtract from the coalition value the values already entered
        val <- v[k] - sum(contributions[i, ], na.rm = TRUE)
        # Save it in the column determined by the corresponding order
        contributions[i, actualOrder[j]]  <- val
        
      }
      
    }
    
    # Computes average for each player
    shapleyVal <- colSums(contributions) / nrow(orders)
    
    # Returns contributions matrix and Shapley value
    return(list(contributions = contributions, value = shapleyVal))
    
  }
  
}
#-----------------------------------------------------------------------------#