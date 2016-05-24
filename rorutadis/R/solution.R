#### SOLVING MODEL

#' @import Rglpk
extremizeVariable <- function(constraints, variableIndex, maximize) {
  obj <- rep(0, ncol(constraints$lhs))
  obj[variableIndex] <- 1  
  Rglpk_solve_LP(obj, constraints$lhs, constraints$dir, constraints$rhs, max = maximize,
                 types = constraints$types)
}

maximizeEpsilon <- function(model) {
  stopifnot(!is.null(model$epsilonIndex))
  
  return (extremizeVariable(model$constraints, model$epsilonIndex, TRUE))
}

isModelConsistent <- function(model) {            
  ret <- maximizeEpsilon(model)
  
  return (ret$status == 0 && ret$optimum >= RORUTADIS_MINEPS)
}

getThresholdsFromF <- function(model, values) {
  return (values[model$firstThresholdIndex:(model$firstThresholdIndex + model$nrClasses - 2)])
}

getAssignmentsFromF <- function(model, values, thresholds) {
  nrVariables <- ncol(model$constraints$lhs)
  assignments <- c()
  
  for (i in seq_len(nrow(model$perfToModelVariables))) {
    assignments[i] <- 1
    
    for (h in seq_len(model$nrClasses - 1)) {
      if (sum(ua(i, nrVariables, model$perfToModelVariables) * values) < thresholds[h]) { #todo: consider epsilon
        break
      }
      
      assignments[i] <- assignments[i] + 1
    }
  }
  
  return (assignments)
}

toSolution <- function(model, values) {
  nrVariables <- ncol(model$constraints$lhs)
  nrAlternatives <- nrow(model$perfToModelVariables)
  nrCriteria <- ncol(model$perfToModelVariables)
    
  stopifnot(length(values) == nrVariables)
  
  # thresolds
  
  thresholds <- getThresholdsFromF(model, values)
  
  # assignments
  
  assignments <- getAssignmentsFromF(model, values, thresholds)
  
  # epsilon
  
  epsilon <- NULL
  
  if (!is.null(model$epsilonIndex)) {
    epsilon <- values[model$epsilonIndex]
  } else {
    epsilon <- RORUTADIS_MINEPS
  }
  
  # vf
  
  vf <- list()
  
  for (j in seq_len(nrCriteria)) {
    nrValues <- length(model$criterionValues[[j]])
    
    if (model$generalVF[j]) {
      x <- sapply(model$criterionValues[[j]], function(w) { w$value })
    } else {
      firstValue <- model$criterionValues[[j]][[1]]$value
      lastValue <- model$criterionValues[[j]][[length(model$criterionValues[[j]])]]$value
      intervalLength <- (lastValue - firstValue) / (model$chPoints[j] - 1)
      
      x <- c(firstValue,
             unlist(sapply(seq_len(model$chPoints[j] - 2), function(w) { firstValue + intervalLength * w })),
             lastValue)
    }
    
    y <- values[model$firstChPointVariableIndex[j] : (model$firstChPointVariableIndex[j] + model$chPoints[j] - 2)]
        
    if (model$criterionPreferenceDirection[j] == "g") {
      y <- c(0, y)
    } else {
      y <- c(y, 0)
    }
    
    vf[[j]] <- cbind(x, y)
  }
  
  # alternative values
  alternativeValues <- matrix(nrow=nrAlternatives, ncol=nrCriteria)
  
  for (i in seq_len(nrAlternatives)) {
    for (j in seq_len(nrCriteria)) {
      alternativeValues[i, j] <- 0
      
      for (k in seq_len(length(model$perfToModelVariables[[i, j]]))) {
        alternativeValues[i, j] <- alternativeValues[i, j] + values[model$perfToModelVariables[[i, j]][[k]][1]] * model$perfToModelVariables[[i, j]][[k]][2]
      }
    }
  }
  
  return (list(
    vf = vf,
    thresholds = thresholds,
    assignments = assignments,
    alternativeValues = alternativeValues,
    solution = values,
    epsilon = epsilon,
    generalVF = model$generalVF
    ))
}

#### SOLUTION

#' Get thresholds
#'
#' This function extracts values of thresholds from solution.
#' 
#' Function is deprecated. Solution already contains thresholds.
#' 
#' @param problem Problem whose model was solved.
#' @param solution Result of model solving (e.g. result of
#' \code{\link{findRepresentativeFunction}} or \code{\link{investigateUtility}}).
#' @return Vector containing \code{h-1} thresholds from \code{t_1} to \code{t_h-1}
#' where \code{t_p-1} is lower threshold of class \code{C_p} and \code{h} is
#' number of classes.
#' @export
getThresholds <- function(problem, solution) {
  .Deprecated()
  return (solution$thresholds)
}

#' Get marginal utilities
#'
#' This function extracts alternatives marginal values from model solution.
#' 
#' Function is deprecated. Solution already contains marginal utilities.
#' 
#' @param problem Problem whose model was solved.
#' @param solution Result of model solving (e.g. result of
#' \code{\link{findRepresentativeFunction}} or \code{\link{investigateUtility}}).
#' @return A \emph{n} x \emph{m} matrix containing marginal values of \code{n} alternatives
#' on \code{m} criteria.
#' @export
getMarginalUtilities <- function(problem, solution) {
  .Deprecated()
  return (solution$alternativeValues)
}

#' Get characteristic points
#'
#' This function extracts values of characteristic points from model solution.
#' 
#' Function is deprecated. Solution already contains characteristic points.
#' 
#' @param problem Problem whose model was solved.
#' @param solution Result of model solving (e.g. result of
#' \code{\link{findRepresentativeFunction}} or \code{\link{investigateUtility}}).
#' @return List of \code{m} matrices for each of \code{m} criteria.
#' Each row \code{c(g, u)} of each matrix contains coordinates of a single
#' characteristic point, where \code{g} - evaluation on corresponding criterion,
#' \code{u} - marginal utility.
#' @export
getCharacteristicPoints <- function(problem, solution) {
  .Deprecated()
  return (solution$vf)
}


#' Get assignments
#'
#' This function returns assignments for given model solution.
#' 
#' Function is deprecated. Solution already contains assignments.
#' 
#' @param problem Problem whose model was solved.
#' @param solution Result of model solving (e.g. result of
#' \code{\link{findRepresentativeFunction}} or \code{\link{investigateUtility}}).
#' @return Vector of alternative assignments. Each element contains an index
#' of a class that corresponding alternative was assigned to.
#' @export
getAssignments <- function(problem, solution) {
  .Deprecated()
  return (solution$assignments)
}

