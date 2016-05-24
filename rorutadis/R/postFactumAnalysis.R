#### HELPERS


isModelConsistentForRho <- function(alternative, atLeastToClass, criteria, necessary,
                                    improvement, limits, problem, rho) {
  original <- problem$perf[alternative, ]
  
  if (improvement) {
    for (criterion in 1:ncol(problem$perf)) {
      if (criteria[criterion]) {
        problem$perf[alternative, criterion] <- min(original[criterion] * rho,
                                                    limits[criterion])
      }
    }
  } else {
    # deterioration
    for (criterion in 1:ncol(problem$perf)) {
      if (criteria[criterion]) {
        problem$perf[alternative, criterion] <- max(original[criterion] * rho,
                                                    limits[criterion])
      }
    }
  }
  
  model <- buildModel(problem, T)
  
  if (atLeastToClass > 1) {
    if (necessary) {
      model$constraints <- combineConstraints(model$constraints,
                                              buildUBAssignmentsConstraint(alternative, atLeastToClass - 1, model))
    } else {
      model$constraints <- combineConstraints(model$constraints,
                                              buildLBAssignmentsConstraint(alternative, atLeastToClass, model))
    }
  }
  
  return (isModelConsistent(model))
}


#### POST FACTUM ANALYSIS

#' Post factum analysis: deteriorate assignment
#'
#' This function checks how much an alternative evaluations can be deteriorated
#' so that that alternative would stay possibly (or necessarily)
#' in at least some specific class. Deterioration is based on minimization value of
#' \code{rho} in multiplication of an alternative evaluations on selected
#' criteria by value \code{rho} (where \code{0 < rho <= 1}). 
#' \strong{Note!} This function works for problems with only non-negative
#' alternative evaluations.
#' @param alternative An alternative for assignment deterioration.
#' @param atLeastToClass An assignment to investigate.
#' @param criteriaManipulability Vector containing a logical value for each criterion.
#' Each value denotes whether multiplying by \code{rho} on corresponding criterion is allowed or not.
#' At least one criterion has to be available for that manipulation.
#' @param necessary Whether necessary or possible assignment is considered.
#' @param problem Problem for which deterioration will be performed.
#' @return Value of \code{rho} or \code{NULL} if given assignment is not possible
#' in any scenario.
#' @seealso
#' \code{\link{improveAssignment}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.5), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' rho <- deteriorateAssignment(4, 1, c(TRUE, TRUE), FALSE, problem)
#' @export
deteriorateAssignment <- function(alternative,
                                  atLeastToClass,
                                  criteriaManipulability,
                                  necessary,
                                  problem) {
  if (min(problem$perf) < 0) {
    stop("Function deteriorateAssignment works only for cases with non-negative evaluations.")
  }
  
  if (length(which(criteriaManipulability)) == 0) {
    stop("At least one criterion has to be active for manipulation.")
  }
  
  if (length(which(problem$criteria == 'c')) != 0) {
    stop("Function deteriorateAssignment works only for cases with 'gain' criteria.")
  }
  
  if (sum(problem$perf[alternative, which(criteriaManipulability)]) == 0) {
    stop("All performances of the alternative on selected criteria are equal to zero. Analysis cannot be performed.")
  }
  
  if (necessary) {
    if (is.null(deteriorateAssignment(alternative, atLeastToClass,
                                      criteriaManipulability, FALSE, problem)))
      return (NULL)
  }
  
  limits <- c()
  original <- problem$perf[alternative, ]
  
  for (criterion in 1:ncol(problem$perf)) {
    if (criteriaManipulability[criterion]) {
      limits <- c(limits, min(problem$perf[, criterion]))
    }
    else {
      limits <- c(limits, problem$perf[alternative, criterion])
    }
  }
  
  minimums <- limits
  
  if (any(criteriaManipulability == FALSE)) {
    minimums <- minimums[-which(criteriaManipulability == FALSE)]
    original <- original[-which(criteriaManipulability == FALSE)]
  }
  
  zeroZero <- which(minimums == 0 & original == 0)
  
  if (length(zeroZero) > 0) {
    minimums <- minimums[-zeroZero]
    original <- original[-zeroZero]
  }
  
  left <- 1
  if (length(original) > 0) {
    left <- min(minimums / original)
  }
  right <- 1
  current <- 1# not left + (right - left)/ 2 due to an immediate and accurate response if deterioration is not possible
  valueForLastTrue <- NULL
  
  repeat {
    if (isModelConsistentForRho(alternative, atLeastToClass, criteriaManipulability, necessary,
                                FALSE, limits, problem, current)) {
      if (RORUTADIS_VERBOSE) print (paste("Model is consistent for rho =", current))
      
      if (necessary) {
        if (is.null(valueForLastTrue)) {
          valueForLastTrue <- current
        }
        else if (valueForLastTrue < current) {
          valueForLastTrue <- current
        }
        
        left <- current
        current <- current + (right - current) / 2
      }
      else {
        if (is.null(valueForLastTrue)) {
          valueForLastTrue <- current
        }
        else if (valueForLastTrue > current) {
          valueForLastTrue <- current
        }
        
        right <- current
        current <- left + (current - left) / 2
      }
    }
    else {
      if (RORUTADIS_VERBOSE) print (paste("Model is inconsistent for rho =", current))
      
      if (necessary) {
        right <- current
        current <- left + (current - left) / 2
      }
      else {
        left <- current
        current <- current + (right - current) / 2
      }
    }
    
    if (right - left < RORUTADIS_MINEPS) {
      break;
    }
  }
  
  return (valueForLastTrue)
}


#' Post factum analysis: improve assignment
#' 
#' This function calculates minimal \code{rho} by which
#' alternative evaluations on selected criteria have to be multiplied for that alternative to be 
#' possibly (or necessarily) assigned to at least some specific class
#' (\code{rho >= 1}). 
#' \strong{Note!} This function works for problems with only non-negative
#' alternative evaluations.
#' @param alternative An alternative for assignment improvement.
#' @param atLeastToClass Desired assignment.
#' @param criteriaManipulability Vector containing a logical value for each criterion.
#' Each value denotes whether multiplying by \code{rho} on corresponding criterion is allowed or not.
#' At least one criterion has to be available for that manipulation.
#' @param necessary Whether necessary or possible assignment is considered.
#' @param problem Problem for which improvement will be performed.
#' @return Value of \code{rho} or \code{NULL} if given assignment is not possible
#' in any scenario.
#' @seealso
#' \code{\link{deteriorateAssignment}}
#' @examples
#' perf <- matrix(c(8, 2, 1, 7, 0.5, 0.9, 0.4, 0.5), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsUB(problem, c(1, 2), c(2, 3))
#' 
#' # a_1 dominates a_4 and a_1 is assigned at most to class C_2
#' # How many times evaluations of a_4 should be improved
#' # that a_4 will be assigned possibly to class C_3?
#' rho <- improveAssignment(4, 3, c(TRUE, TRUE), FALSE, problem)
#' @export
improveAssignment <- function(alternative,
                              atLeastToClass,
                              criteriaManipulability,
                              necessary,
                              problem) {
  if (min(problem$perf) < 0) {
    stop("Function improveAssignment works only for cases with non-negative evaluations.")
  }
  
  if (length(which(criteriaManipulability)) == 0) {
    stop("At least one criterion has to be active for manipulation.")
  }
  
  if (length(which(problem$criteria == 'c')) != 0) {
    stop("Function improveAssignment works only for cases with 'gain' criteria.")
  }
  
  if (sum(problem$perf[alternative, which(criteriaManipulability)]) == 0) {
    stop("All performances of the alternative on selected criteria are equal to zero. Analysis cannot be performed.")
  }
      
  if (necessary) {
    if (is.null(improveAssignment(alternative, atLeastToClass, criteriaManipulability, FALSE, problem)))
      return (NULL)
  }
  
  limits <- c()
  original <- problem$perf[alternative, ]
  
  for (criterion in 1:ncol(problem$perf)) {
    if (criteriaManipulability[criterion]) {
      limits <- c(limits, max(problem$perf[, criterion]))
    } else {
      limits <- c(limits, problem$perf[alternative, criterion])
    }
  }
  
  maximums <- limits
  
  if (any(criteriaManipulability == FALSE)) {
    maximums <- maximums[-which(criteriaManipulability == FALSE)]
    original <- original[-which(criteriaManipulability == FALSE)]
  }
  
  zero <- which(original == 0)
  
  if (length(zero) > 0) {
    maximums <- maximums[-zero]
    original <- original[-zero]
  }
  
  left <- 1
  right <- 1
  if (length(original) > 0) {
    right <- max(maximums / original)
  }  
  current <- 1# not left + (right - left)/ 2 due to an immediate and accurate response if improvement is not possible
  valueForLastTrue <- NULL
  
  repeat {
    if (isModelConsistentForRho(alternative, atLeastToClass, criteriaManipulability, necessary,
                                TRUE, limits, problem, current)) {
      if (RORUTADIS_VERBOSE) print (paste("Model is consistent for rho =", current))
      
      if (necessary) {
        if (is.null(valueForLastTrue)) {
          valueForLastTrue <- current
        }
        else if (valueForLastTrue < current) {
          valueForLastTrue <- current
        }
        
        left <- current
        current <- current + (right - current) / 2
      }
      else {
        if (is.null(valueForLastTrue)) {
          valueForLastTrue <- current
        }
        else if (valueForLastTrue > current) {
          valueForLastTrue <- current
        }
        
        right <- current
        current <- left + (current - left) / 2
      }
    }
    else {
      if (RORUTADIS_VERBOSE) print (paste("Model is inconsistent for rho =", current))
      
      if (necessary) {
        right <- current
        current <- left + (current - left) / 2
      }
      else {
        left <- current
        current <- current + (right - current) / 2
      }
    }
    
    if (right - left < RORUTADIS_MINEPS) {
      break;
    }
  }
  
  return (valueForLastTrue)
}

#' Post factum analysis: check how much utility is missing
#'
#' This function calculates missing value of an alternative utility for
#' that alternative to be possibly (or necessarily) assigned to at least some
#' specific class.
#' 
#' @param alternative An alternative index.
#' @param atLeastToClass An assignment to investigate.
#' @param necessary Whether necessary or possible assignment is considered.
#' @param problem Problem for investigation.
#' @return List with named elements:
#' \itemize{
#' \item \code{ux} - value of missing utility,
#' \item \code{solution} - result of solving model. It can be used for further
#' computations (\code{\link{getAssignments}}, \code{\link{getThresholds}}, \code{\link{getMarginalUtilities}},
#' \code{\link{getCharacteristicPoints}}).
#' }
#' \code{NULL} is returned if given assignment is not possible.
#' @seealso
#' \code{\link{getMarginalUtilities}}
#' \code{\link{getCharacteristicPoints}}
#' \code{\link{getThresholds}}
#' \code{\link{improveAssignment}}
#' @examples
#' perf <- matrix(c(8, 2, 1, 7, 0.5, 0.9, 0.4, 0.5), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsUB(problem, c(1, 2), c(2, 3))
#' 
#' result <- investigateUtility(4, 3, FALSE, problem)
#' @export
investigateUtility <- function(alternative, atLeastToClass, necessary, problem) {
  stopifnot(atLeastToClass > 1)
  stopifnot(atLeastToClass <= problem$nrClasses)
  
  model <- buildModel(problem, FALSE)
  model$constraints <- addVarialbesToModel(model$constraints, c("C", "C"))
  uxIndex <- ncol(model$constraints$lhs) - 1
  
  if (atLeastToClass > 1) {
    if (necessary) {
      constr <- buildUBAssignmentsConstraint(alternative, atLeastToClass - 1, model)
      constr$rhs <- 0
    } else {
      constr <- buildLBAssignmentsConstraint(alternative, atLeastToClass, model)
    }
    
    constr$lhs[uxIndex] <- 1
    constr$lhs[uxIndex + 1] <- -1
    
    model$constraints <- combineConstraints(model$constraints, constr)
  }
  
  obj <- rep(0, ncol(model$constraints$lhs))
  obj[uxIndex] <- 1
  obj[uxIndex + 1] <- -1
  ret <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                        max = necessary, types = model$constraints$types)

  if (ret$status == 0) {
    ux <- 0
    
    if (ret$optimum > 0) {
      ux <- ret$optimum
    }
    
    return (list(ux = ux, solution = toSolution(model, ret$solution)))
  } else {
    return (NULL)
  }
}
