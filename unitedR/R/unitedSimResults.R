#' @include summary.R
#' @include unitedSimClass.R
NULL

###############################################
# --------------------------------------------#
# Class unitedSimResults                      #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the unitedSim class
# 
# @param object \code{S4} object of the \code{unitedSim} class
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateUnitedSimResults <- function(object) {
  errors <- character()
  
  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for unitedSimResults
# --------------------------------------------

# unitedSim class
setClass("unitedSimResults",
         slots = c(games = "list"),
         validity = validateUnitedSimResults)


# --------------------------------------------
# Methods for unitedSimResults
# --------------------------------------------

#' @rdname summary
setMethod("summary", signature(object = "unitedSimResults"), function(object) {
    summaries<- lapply(object@games, function(l) {
        summary(l)
      }
    )  
    summaries <- do.call("rbind", summaries)
    rownames(summaries) <- 1:length(object@games)
    summaries
  }
)

# --------------------------------------------
# Show function for unitedSimResults
# --------------------------------------------

setMethod("show", "unitedSimResults", function(object) {
    # iterate through all slots of the object
    cat("\n")
    cat("The used lineup home\n")
    lineupHome <- toString(getLineup(object@games[[1]]@home))
    lineupHome <- gsub(", ", "-", lineupHome)
    cat("\t", lineupHome)
    cat("\nwas compared to the following away lineups\n")
    D <- summary(object)
    D$home <- NULL
    D$averageTrainingPointsAway <- NULL
    D$averagePointsAway <- NULL
    D$winProbabilityAway <- NULL
    print(D, row.names = FALSE)
    
    #lineupsAway <- sapply(object@games, function(l) {
    #    lineup <- toString(getLineup(l@away))
    #    lineup <- gsub(", ", "-", lineup)
    #  }
    #)  
    #winProbability <- sapply(object@games, function(l) {
    #    round(l@winProbabilityHome, digits = 4)
    #  }
    #)  
    #for (i in 1:length(lineupsAway)) {
    #  cat("\t", lineupsAway[i], "\twinProbability:", winProbability[i])
    #  cat("\n")
    #}
  }  
)