#' @include formation.R
NULL

###############################################
# --------------------------------------------#
# penaltyGoalsProb                            #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Main Function for computing goals by penalties
# --------------------------------------------

#' Computing goals by united
#' 
#' Computes the distribution function of possible goals by penalties.
#' 
#' @inheritParams overview
#'
#' @return 
#' A \code{data.frame} with two columns: the possible goals and the probability
#' for achieving this number of goals.
#'
#' @export
penaltyGoalsProb <- function(posPenalties, penaltyGoalProb) {
  stopifnot(posPenalties >= 0, posPenalties < 12, round(posPenalties) == posPenalties)
  stopifnot(penaltyGoalProb >= 0, penaltyGoalProb <= 1)
  # vector for the probabilities of the goals which can be scored by penalties
  probs <- numeric(posPenalties + 1)
  # i variable for scored goals
    for (i in 0:posPenalties) {
      prob <- 0
      # j variable of achieved penalties
      for (j in i:posPenalties) {
        probGoals <- choose(j, i) * penaltyGoalProb^i * (1-penaltyGoalProb)^(j-i)
        prob <- prob + choose(posPenalties, j) * 0.1^j * 0.9^(posPenalties-j) * probGoals
      }
      probs[i+1] <- prob
    }
  data.frame(goals = 0:posPenalties, probability = probs)
}
  
