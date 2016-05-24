###############################################
# --------------------------------------------#
# Class unitedSim                             #
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
validateUnitedSim <- function(object) {
  errors <- character()
  
  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for unitedSim
# --------------------------------------------

# unitedSim class
setClass("unitedSim",
         slots = c(results = "data.frame", averageTrainingPointsHome = "numeric",
                   averageTrainingPointsAway = "numeric", averagePointsHome = "numeric",
                   averagePointsAway = "numeric", winProbabilityHome = "numeric",
                   winProbabilityAway = "numeric", tiedProbability = "numeric",
                   home = "formation", away = "formation"),
         validity = validateUnitedSim)


# --------------------------------------------
# Class definition for unitedSimR
# --------------------------------------------

# Randomization paramters generic
setClass("unitedSimR",
         slots = c(r = "numeric"),
         contains = "unitedSim")

# --------------------------------------------
# Methods for unitedSim
# --------------------------------------------

#' @rdname summary
setMethod("summary", signature(object = "unitedSim"), function(object) {
    names <- slotNames(object)
    names <- names[!(names %in% c("results", "home", "away"))]
    value <- numeric(length(names) + 2)
    lineupHome <- toString(getLineup(object@home))
    value[1] <- lineupHome <- gsub(", ", "-", lineupHome)
    lineupAway <- toString(getLineup(object@away))
    value[2] <- gsub(", ", "-", lineupAway)
    for (i in 1:length(names)) {
      value[i+2] <- slot(object, names[i])
    }
    # paste a data.frame
    D <- as.data.frame(t(data.frame(value)))
    colnames(D) <- c("home", "away", names)
    rownames(D) <- 1
    for (i in 1:2) {
      D[ ,i] <- as.character(D[ ,i])
    }
    for (i in 3:ncol(D)) {
      D[ ,i] <- as.numeric(as.character(D[ ,i]))
    }
    D
  }
)


# --------------------------------------------
# Show function for unitedSim
# --------------------------------------------

setMethod("show", "unitedSim", function(object) {
    # iterate through all slots of the object
    cat("\n")
    cat("Used lineup home was:\n")
    lineupHome <- toString(getLineup(object@home))
    lineupHome <- gsub(", ", "-", lineupHome)
    cat("\t", lineupHome)
    cat("\n")
    cat("Used lineup away was:\n")
    lineupAway <- toString(getLineup(object@away))
    lineupAway <- gsub(", ", "-", lineupAway)
    cat("\t", lineupAway)
    cat("\n\nThe key statistcs are:\n")
    names <- slotNames(object)
    names <- names[!(names %in% c("results", "home", "away"))]
    for (name in names) {
      cat("\t", name, "=", slot(object, name), "\n")
    }
    cat("\nThe most probable results are:\n") 
    if (nrow(object@results) >= 6) {
      print(round(object@results[1:6, 1:4], digits = 3), row.names = FALSE)
    } else {
      print(round(object@results[, 1:4], digits = 3), row.names = FALSE)
    }  
    cat("\n")
  }  
)


# --------------------------------------------
# Show function for unitedSim
# --------------------------------------------

setMethod("show", "unitedSimR", function(object) {
    # iterate through all slots of the object
    cat("\n")
    cat("Used lineup home was:\n")
    lineupHome <- toString(getLineup(object@home))
    lineupHome <- gsub(", ", "-", lineupHome)
    cat("\t", lineupHome)
    cat("\n")
    cat("Used lineup away was:\n")
    lineupAway <- toString(getLineup(object@away))
    lineupAway <- gsub(", ", "-", lineupAway)
    cat("\t", lineupAway)
    cat("\n\nThe key statistcs based on", object@r,"simulations are:\n")
    names <- slotNames(object)
    names <- names[!(names %in% c("results", "home", "away", "r"))]
    for (name in names) {
      cat("\t", name, "=", slot(object, name), "\n")
    }
    cat("\nThe most probable results are:\n") 
    if (nrow(object@results) >= 6) {
      print(round(object@results[1:6, 1:4], digits = 3), row.names = FALSE)
    } else {
      print(round(object@results[, 1:4], digits = 3), row.names = FALSE)
    }  
    cat("\n")
  }  
)
