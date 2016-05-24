#' @include formation.R
#' @include penaltyGoalsProb.R
#' @include unitedSimClass.R
NULL

utils::globalVariables(c("goalsHome", "goalsAway", "probability"))

###############################################
# --------------------------------------------#
# unitedSimOne                                #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Main Function for simulation line ups
# --------------------------------------------

#' Simulating a formation
#' 
#' Simulates a formation against another formation.
#' 
#' @inheritParams overview
#'
#' @return Creates an object of the \code{unitedSim} class.
#' 
#' @examples 
#' home <- formation(10, NA, c(7,5,3), c(8,8), c(10,10,8))
#' away <- formation(5, 8, c(8,8), c(10,10), c(10,10,10), 
#'  hardness = c(0,0,0,0,1))
#' set.seed(123)
#' unitedSimOne(home, away)
#' # you can even simulated the game
#' unitedSimOne(home, away, r = 100)
#' 
#' @export
unitedSimOne <- function(home, away, r, preventGoalGK = 1/14, preventGoalSW = 1/15) {
  stopifnot(validObject(home), validObject(away), is(home, "formation"), 
            is(home, "formation"), is.numeric(preventGoalGK), is.numeric(preventGoalSW))
  if (preventGoalGK >= 1/13) stop("preventGoalGK must be smaller than 1/13.")
  if (preventGoalGK < 0) stop("preventGoalGK must be greater than zero.")
  if (preventGoalSW >= 1/13) stop("preventGoalSW must be smaller than 1/13.")
  if (preventGoalSW < 0) stop("preventGoalSW must be greater than zero.")
  if (missing(r)) {
    if (sum(home@hardness) > 1 || sum(away@hardness > 1)) {
      warning("It is recommended to simulate hardness and penalties, calculations are exact for one possible lineup.")
    }
    
    homeLineup <- getLineup(home)
    awayLineup <- getLineup(away)
  
    # simulate red cards
    homeLineupSim <- simRedCard(home, homeLineup)
    awayLineupSim <- simRedCard(away, awayLineup)
  
    chancesHome <- round((homeLineupSim[3:5] - awayLineupSim[5:3] - 
                          c(0, 0, awayLineupSim[2])) * c(1/4, 1/2, 1) + 0.00001)
    chancesHome <- sum(chancesHome[chancesHome > 0])
  
    chancesAway <- round((awayLineupSim[3:5] - homeLineupSim[5:3] - 
                          c(0, 0, homeLineupSim[2])) * c(1/4, 1/2, 1) + 0.00001)
    chancesAway <- sum(chancesAway[chancesAway > 0])
  
    # possible penalties home
    posPenaltiesHome <- sum(away@hardness)
    posPenaltiesAway <- sum(home@hardness)
  
    # probability of a goal of a penalty
    penaltyProbGoalHome <- 1 - (awayLineupSim[1] * 0.05)
    penaltyProbGoalAway <- 1 - (homeLineupSim[1] * 0.05)
  
    # probability distribution of all possible goals by penalties for both teams
    goalsPenaltyDistrHome <- penaltyGoalsProb(posPenaltiesHome, penaltyProbGoalHome)
    goalsPenaltyDistrAway <- penaltyGoalsProb(posPenaltiesAway, penaltyProbGoalAway)
  
    # possible allocations of the penalties
    penaltyAllocations <- expand.grid(home = 0:posPenaltiesHome, away = 0:posPenaltiesAway)
    penaltyAllocations$probability <- apply(penaltyAllocations, 1, function(x) {
        probHome <- goalsPenaltyDistrHome$probability[which(goalsPenaltyDistrHome$goals == x[1])]
        probAway <- goalsPenaltyDistrAway$probability[which(goalsPenaltyDistrAway$goals == x[2])]
        probHome * probAway
      }
    )
  
    # compute all possible results
    allResults <- expand.grid(goalsHome = 0:chancesHome, goalsAway = 0:chancesAway)
    # compute probability of a goal
    probGoalAway <- (1-homeLineupSim[1] * preventGoalGK) * (1-homeLineupSim[2] * preventGoalSW)
    probGoalHome <- (1-awayLineupSim[1] * preventGoalGK) * (1-awayLineupSim[2] * preventGoalSW)
    # compute probability of results
    allResults$probability <- apply(allResults, 1, function(x) {
        probHome <- dbinom(x[1], chancesHome, prob = probGoalHome)
        probAway <- dbinom(x[2], chancesAway, prob = probGoalAway)
        probHome * probAway
      }
    )
  
    # list for all possible results including penalties
    allResWithPen <- vector("list", nrow(penaltyAllocations))
    for (i in 1:nrow(penaltyAllocations)) {
      resultsPenalty <- allResults
      resultsPenalty$goalsHome <- resultsPenalty$goalsHome + penaltyAllocations$home[i]
      resultsPenalty$goalsAway <- resultsPenalty$goalsAway + penaltyAllocations$away[i]
      resultsPenalty$probability <- resultsPenalty$probability * penaltyAllocations$probability[i]
      allResWithPen[[i]] <- resultsPenalty
    }
  
    # merge the data.frames
    finalPossibleResults <- do.call("rbind", allResWithPen)
    
    finalPossibleResults <- ddply(finalPossibleResults, .(goalsHome, goalsAway), summarize, 
                                probability = sum(probability))
  
  
    # sort results by probability of apprearance
    finalPossibleResults <- finalPossibleResults[order(finalPossibleResults$probability, decreasing = TRUE), ]
    # add a cumsum of the probabilities
    finalPossibleResults$cumsumProb <- cumsum(finalPossibleResults$probability)
    # add points for home
    finalPossibleResults$pointsHome <- ifelse(finalPossibleResults$goalsHome > finalPossibleResults$goalsAway, 3, 
                                            ifelse(finalPossibleResults$goalsHome == finalPossibleResults$goalsAway, 1,
                                                   0))
    # add points for away
    finalPossibleResults$pointsAway <- ifelse(finalPossibleResults$goalsHome < finalPossibleResults$goalsAway, 3, 
                                            ifelse(finalPossibleResults$goalsHome == finalPossibleResults$goalsAway, 1, 
                                                   0))
    # add training points (TP) for home
    finalPossibleResults$tpHome <- ifelse(finalPossibleResults$goalsHome > finalPossibleResults$goalsAway, 1, 
                                        ifelse(finalPossibleResults$goalsHome == finalPossibleResults$goalsAway, 0.5, 
                                               0))
    # add traings points (TP) for away
    finalPossibleResults$tpAway <- ifelse(finalPossibleResults$goalsHome < finalPossibleResults$goalsAway, 1, 
                                        ifelse(finalPossibleResults$goalsHome == finalPossibleResults$goalsAway, 0.5, 
                                               0))
    # output
    output <- new("unitedSim", 
                results = finalPossibleResults, 
                averageTrainingPointsHome = round(sum(finalPossibleResults$tpHome * finalPossibleResults$probability), digits = 4),
                averageTrainingPointsAway = round(sum(finalPossibleResults$tpAway * finalPossibleResults$probability), digits = 4), 
                averagePointsHome = round(sum(finalPossibleResults$pointsHome * finalPossibleResults$probability), digits = 4),
                averagePointsAway = round(sum(finalPossibleResults$pointsAway * finalPossibleResults$probability), digits = 4),
                winProbabilityHome = round(sum((finalPossibleResults$pointsHome == 3) * finalPossibleResults$probability), digits = 4),
                winProbabilityAway = round(sum((finalPossibleResults$pointsAway == 3) * finalPossibleResults$probability), digits = 4), 
                tiedProbability = round(sum((finalPossibleResults$pointsAway == 1) * finalPossibleResults$probability), digits = 4),
                home = home,
                away = away)
  
    return(output)
    # include red Cards when simulating
  } else {
    stopifnot(is.numeric(r), round(r) == r, length(r) == 1)
    homeLineup <- getLineup(home)
    awayLineup <- getLineup(away)
    simulatedResults <- t(sapply(1:r, function(x) { 
              # simulate red cards
              homeLineupSim <- simRedCard(home, homeLineup)
              awayLineupSim <- simRedCard(away, awayLineup)
              
              chancesHome <- round((homeLineupSim[3:5] - awayLineupSim[5:3] - 
                                      c(0, 0, awayLineupSim[2])) * c(1/4, 1/2, 1))
              chancesHome <- sum(chancesHome[chancesHome > 0])
              
              chancesAway <- round((awayLineupSim[3:5] - homeLineupSim[5:3] - 
                                      c(0, 0, homeLineupSim[2])) * c(1/4, 1/2, 1))
              chancesAway <- sum(chancesAway[chancesAway > 0])
              
              #  penalties home
              penaltiesHome <- rbinom(1, sum(away@hardness), 0.1)
              penaltiesAway <- rbinom(1, sum(home@hardness), 0.1)
              
              # probability of a goal by penalty
              penaltyProbGoalHome <- 1 - (awayLineupSim[1] * 0.05)
              penaltyProbGoalAway <- 1 - (homeLineupSim[1] * 0.05)
              
              # probability of a goal with a chance
              probGoalAway <- (1-homeLineupSim[1] * preventGoalGK) * (1-homeLineupSim[2] * preventGoalSW)
              probGoalHome <- (1-awayLineupSim[1] * preventGoalGK) * (1-awayLineupSim[2] * preventGoalSW)
              
              # simulate the game
              penaltyGoalsHome <- rbinom(1, penaltiesHome, penaltyProbGoalHome)
              penaltyGoalsAway <- rbinom(1, penaltiesAway, penaltyProbGoalAway)
              goalsHomeGame <- rbinom(1, chancesHome, probGoalHome)
              goalsAwayGame <- rbinom(1, chancesAway, probGoalAway)
              c(penaltyGoalsHome + goalsHomeGame, penaltyGoalsAway + goalsAwayGame)
      }
    ))
    simulatedResults <- as.data.frame(simulatedResults)
    colnames(simulatedResults) <- c("goalsHome", "goalsAway")
    simulatedResults$probability <- 1/r
    
    simulatedResults <- ddply(simulatedResults, .(goalsHome, goalsAway), summarize, 
                                  probability = sum(probability))
    
    # sort results by probability of apprearance
    simulatedResults <- simulatedResults[order(simulatedResults$probability, decreasing = TRUE), ]
    # add a cumsum of the probabilities
    simulatedResults$cumsumProb <- cumsum(simulatedResults$probability)
    # add points for home
    simulatedResults$pointsHome <- ifelse(simulatedResults$goalsHome > simulatedResults$goalsAway, 3, 
                                              ifelse(simulatedResults$goalsHome == simulatedResults$goalsAway, 1,
                                                     0))
    # add points for away
    simulatedResults$pointsAway <- ifelse(simulatedResults$goalsHome < simulatedResults$goalsAway, 3, 
                                              ifelse(simulatedResults$goalsHome == simulatedResults$goalsAway, 1, 
                                                     0))
    # add training points (TP) for home
    simulatedResults$tpHome <- ifelse(simulatedResults$goalsHome > simulatedResults$goalsAway, 1, 
                                          ifelse(simulatedResults$goalsHome == simulatedResults$goalsAway, 0.5, 
                                                 0))
    # add traings points (TP) for away
    simulatedResults$tpAway <- ifelse(simulatedResults$goalsHome < simulatedResults$goalsAway, 1, 
                                          ifelse(simulatedResults$goalsHome == simulatedResults$goalsAway, 0.5, 
                                                 0))
    # output
    output <- new("unitedSimR", 
                  results = simulatedResults, 
                  averageTrainingPointsHome = round(sum(simulatedResults$tpHome * simulatedResults$probability), digits = 4),
                  averageTrainingPointsAway = round(sum(simulatedResults$tpAway * simulatedResults$probability), digits = 4), 
                  averagePointsHome = round(sum(simulatedResults$pointsHome * simulatedResults$probability), digits = 4),
                  averagePointsAway = round(sum(simulatedResults$pointsAway * simulatedResults$probability), digits = 4),
                  winProbabilityHome = round(sum((simulatedResults$pointsHome == 3) * simulatedResults$probability), digits = 4),
                  winProbabilityAway = round(sum((simulatedResults$pointsAway == 3) * simulatedResults$probability), digits = 4), 
                  tiedProbability = round(sum((simulatedResults$pointsAway == 1) * simulatedResults$probability), digits = 4),
                  r = r,
                  home = home,
                  away = away)
    
    return(output)
  }
}
  