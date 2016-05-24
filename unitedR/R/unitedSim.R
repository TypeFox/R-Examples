#' @include unitedSimOne.R
#' @include unitedSimResults.R
NULL

###############################################
# --------------------------------------------#
# unitedSim                                   #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Main Function for simulation line ups
# --------------------------------------------

#' Simulating a formation
#' 
#' Simulates a formation against another formations (several formations of away are possible).
#' 
#' @inheritParams overview
#' @param ... several objects of the class \code{formation}
#'
#' @return Creates an object of the \code{unitedSim} class.
#' 
#' @examples 
#' home <- formation(10, NA, c(7,5,3), c(8,8), c(10,10,8))
#' away <- formation(5, 8, c(8,8), c(10,10), c(10,10,10), 
#'  hardness = c(0,0,0,0,1))
#' set.seed(123)
#' unitedSim(home, away)
#' # can also be simualated
#' unitedSim(home, away, r = 100)
#' # several away lineups
#' unitedSim(home, away, away)
#' # several away lineups simulated
#' unitedSim(home, away, away, r = 100)
#' 
#' @export
unitedSim <- function(home, ..., r, preventGoalGK = 1/14, preventGoalSW = 1/15) {
  stopifnot(validObject(home), is(home, "formation"))
  formations <- list(...)
  if (!all(sapply(formations, function(x)  is(x, "formation"))))
    stop("Not all ... objects of class formation.")
  if (missing(r)) {
    if (length(formations) == 1) {
      return(unitedSimOne(home, formations[[1]], preventGoalGK = preventGoalGK, preventGoalSW = preventGoalSW))
    } else {
      games <- lapply(formations, function(formation) {
          unitedSimOne(home, formation, preventGoalGK = preventGoalGK, preventGoalSW = preventGoalSW)  
        }
      )
    }
    return(new("unitedSimResults", games = games))
  } else {
    stopifnot(is.numeric(r), round(r) == r, length(r) == 1)
      if (length(formations) == 1) {
        return(unitedSimOne(home, formations[[1]], r = r, preventGoalGK = preventGoalGK, preventGoalSW = preventGoalSW))
      } else {
        games <- lapply(formations, function(formation) {
          unitedSimOne(home, formation, r = r, preventGoalGK = preventGoalGK, preventGoalSW = preventGoalSW)  
        }
      )
    }
    return(new("unitedSimResults", games = games))
  }  
}
