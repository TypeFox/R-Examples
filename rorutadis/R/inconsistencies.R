#' Find inconsistencies in preference information
#'
#' This function finds sets of pieces of preference information that make
#' problem inconsistent.
#' 
#' @param problem Problem to investigate. 
#' @return List of ordered by cardinality sets of indices of preference
#' information that makes problem inconsistent. Use \code{\link{getRestrictions}}
#' on sets to find out related preference information.
#' @examples
#' perf <- matrix(c(1, 2, 2, 1), ncol = 2)
#' problem <- buildProblem(perf, 3, TRUE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsUB(problem, c(1, 1))
#' problem <- addAssignmentsLB(problem, c(2, 2))
#' 
#' checkConsistency(problem) # TRUE
#' 
#' problem <- addAssignmentsLB(problem, c(1, 3)) # added inconsistency
#' 
#' checkConsistency(problem) # FALSE
#' 
#' inconsistencies <- findInconsistencies(problem)
#' 
#' setsOfprefInfo <- lapply(inconsistencies,
#'                          function(x) { getRestrictions(problem, x) })
#' @export
findInconsistencies <- function(problem) {
  if (checkConsistency(problem)) {
    return (list())
  }
  
  result <- list()
  
  model <- buildModel(problem, FALSE)
  
  numberOfPrefereranceInformationRecords <- length(model$prefInfoToConstraints)
  
  if (numberOfPrefereranceInformationRecords == 0) {
    stop("There is no preference information in the problem.")
  }
  
  firstPrefInfoVariable <- ncol(model$constraints$lhs) + 1
  model$constraints <- addVarialbesToModel(model$constraints, rep("B", numberOfPrefereranceInformationRecords))
  
  for (i in seq_len(numberOfPrefereranceInformationRecords)) {
    for (j in seq_len(length(model$prefInfoToConstraints[[i]]))) {
      constraintIndex <- model$prefInfoToConstraints[[i]][j]
      
      if (model$constraints$dir[constraintIndex] == ">=") {
        model$constraints$lhs[constraintIndex, firstPrefInfoVariable + i - 1] <- RORUTADIS_BIGM
      } else if (model$constraints$dir[constraintIndex] == "<=") {
        model$constraints$lhs[constraintIndex, firstPrefInfoVariable + i - 1] <- -RORUTADIS_BIGM
      } else {
        stop("Equalities not supported.")
      }
    }
  }
  
  obj <- c(rep(0, firstPrefInfoVariable - 1), rep(1, numberOfPrefereranceInformationRecords))
  nrOfAvailableRecords <- numberOfPrefereranceInformationRecords
  
  while (nrOfAvailableRecords > 0) {
    solution <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                               max=FALSE, types=model$constraints$types)
    
    if (solution$status == 0) {
      resultVector <- solution$solution[(firstPrefInfoVariable):(firstPrefInfoVariable + numberOfPrefereranceInformationRecords)]
      incSetIndices <- which(resultVector == 1)
      
      result[[length(result) + 1]] <- incSetIndices
      
      lhs <- rep(0, firstPrefInfoVariable + numberOfPrefereranceInformationRecords - 1)
      lhs[firstPrefInfoVariable + incSetIndices - 1] <- 1
      
      model$constraints <- combineConstraints(model$constraints, list(lhs = lhs, dir = "==", rhs = 0))
      
      nrOfAvailableRecords <- nrOfAvailableRecords - length(incSetIndices)
    } else {
      break
    }
  }
  
  return (result)
}
