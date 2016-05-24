#### HELPERS

#' Check problem consistency
#'
#' This function allows to check if preference information is consistent.
#' 
#' @param problem Problem to check. 
#' @return \code{TRUE} if a model of a problem is feasible and \code{FALSE}
#' otherwise.
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' isConsistent <- checkConsistency(problem)
#' @export
checkConsistency <- function(problem) {
  model <- buildModel(problem, TRUE)
  return (isModelConsistent(model))
}

checkRelation <- function(model, alternative, class, necessary) {
  stopifnot(!is.null(model$epsilonIndex))
  
  if (necessary) {
    if (class == 1) {
      additionalConstraints <- buildLBAssignmentsConstraint(alternative, 2, model)
    } else if (class == model$nrClasses) {
      additionalConstraints <- buildUBAssignmentsConstraint(alternative, model$nrClasses - 1, model)
    } else if (class > 1 && class < model$nrClasses) {
      model$constraints <- addVarialbesToModel(model$constraints, c("B", "B"))
      nrVariables <- ncol(model$constraints$lhs)
      
      constrLB <- buildLBAssignmentsConstraint(alternative, class + 1, model)
      constrUB <- buildUBAssignmentsConstraint(alternative, class - 1, model)
      constrRel <- list(lhs = rep(0, nrVariables), dir = "==", rhs = 1)
      
      constrLB$lhs[nrVariables - 1] <- RORUTADIS_BIGM
      constrUB$lhs[nrVariables] <- -RORUTADIS_BIGM
      constrRel$lhs[nrVariables - 1] <- 1
      constrRel$lhs[nrVariables] <- 1
      
      additionalConstraints <- combineConstraints(constrLB, constrUB, constrRel)
    }
  } else { # possible
    additionalConstraints <- NULL
    
    if (class > 1) {
      additionalConstraints <- combineConstraints(additionalConstraints,
                                                  buildLBAssignmentsConstraint(alternative, class, model))
    }
    
    if (class < model$nrClasses) {
      additionalConstraints <- combineConstraints(additionalConstraints,
                                                  buildUBAssignmentsConstraint(alternative, class, model))
    }
  }
  
  model$constraints <- combineConstraints(model$constraints, additionalConstraints)
  optimizedEpsilon <- maximizeEpsilon(model)
  
  if (necessary) {
    return (optimizedEpsilon$status != 0 || optimizedEpsilon$optimum < RORUTADIS_MINEPS)
  } else {
    return (optimizedEpsilon$status == 0 && optimizedEpsilon$optimum >= RORUTADIS_MINEPS)
  }
}

calculateExtremeClassCardinality <- function(model, class, maximum) {
  firstAssignmentVariableIndex <- model$firstThresholdIndex + model$nrClasses - 1
  
  obj <- rep(0, ncol(model$constraints$lhs))
  for(alternative in seq_len(nrow(model$perfToModelVariables))) {
    obj[firstAssignmentVariableIndex + (alternative - 1) * model$nrClasses + class - 1]  <- 1
  }
  
  solution <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                             max = maximum, types = model$constraints$types)
  
  stopifnot(solution$status == 0)
  
  return (solution$optimum)
}

# TRUE if alternative is always assigned to class at least as good as class of referenceAlternative
compareAssignment <- function(model, alternative, referenceAlternative) {
  stopifnot(!is.null(model$epsilonIndex))
  
  if (alternative == referenceAlternative) {
    return (TRUE)
  }
  
  for (class in seq_len(model$nrClasses - 1)) {
    newModel <- model
    
    constraintForI <- buildUBAssignmentsConstraint(alternative, class, newModel)
    constraintForJ <- buildLBAssignmentsConstraint(referenceAlternative, class + 1, newModel)
    
    newModel$constraints <- combineConstraints(newModel$constraints, constraintForI, constraintForJ)
    
    if (isModelConsistent(newModel)) {
      return (FALSE)
    }
  }
  
  return (TRUE)
}

#### MAIN PUBLIC FUNCTIONS

#' Calculate assignments
#'
#' This function calculates possible and necessary assignments.
#'
#' @param problem Problem for which assignments will be calculated.
#' @param necessary Whether necessary or possible assignments. 
#' @return \emph{n} x \emph{p} logical matrix, where each row represents one
#' of \emph{n} alternatives and each column represents one of \emph{p} classes.
#' Element \code{[i, h]} is \code{TRUE} if:
#' \itemize{
#' \item for necessary assignments: alternative \code{a_i} is always assigned to
#' class \code{C_h},
#' \item for possible assignments: alternative \code{a_i} can be assigned
#' to class \code{C_h}.
#' }
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' possibleAssignments <- calculateAssignments(problem, FALSE)
#' necessaryAssignments <- calculateAssignments(problem, TRUE)
#' @export
calculateAssignments <- function(problem, necessary) {
  if (!checkConsistency(problem))
    stop("Model infeasible.")
  
  model <- buildModel(problem, T)
  rel <- matrix(nrow = nrow(problem$perf), ncol = problem$nrClasses)
  
  for (alternative in seq_len(nrow(problem$perf))) {
    for(class in seq_len(problem$nrClasses)) {
      rel[alternative, class] <- checkRelation(model, alternative, class, necessary)
      if (RORUTADIS_VERBOSE) print (paste("relation ", alternative, class, "is", rel[alternative, class]))
    }
  }
  
  if (!is.null(rownames(problem$perf)))
    rownames(rel) <- rownames(problem$perf)
  
  return (rel)
}

#' Compare assignments
#'
#' This function compares assignments.
#'
#' @param problem Problem for which assignments will be compared.
#' @param necessary Whether necessary or possible assignments.
#' @return \emph{n} x \emph{n} logical matrix, where \code{n} is a number of
#' alternatives. Cell \code{[i, j]} is \code{TRUE} if \emph{a_i} is assigned to
#' class at least as good as class of \emph{a_j} for all compatible value
#' functions.
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' resultOfComparison <- compareAssignments(problem)
#' @export
compareAssignments <- function(problem, necessary = TRUE) {
  if (necessary == FALSE) {
    stop("Comparing possible assignments not supported.")
  } 
  
  nrAlternatives <- nrow(problem$perf)
  model <- buildModel(problem, T)
  
  result <- matrix(nrow = nrAlternatives, ncol = nrAlternatives)
  
  for (i in seq_len(nrAlternatives)) {
    for (j in seq_len(nrAlternatives)) {
      result[i, j] <- compareAssignment(model, i, j)
            
      if (RORUTADIS_VERBOSE) print (paste("It is ", result[i, j], " that alternative ",
                                          i, " is always in at least as good class as class of alternative ",
                                          j, ".", sep = ""))
    }
  }
  
  return (result)
}

#' Calculate extreme class cardinalities
#'
#' This function calculates minimal and maximal possible cardinality
#' of each class.
#'
#' @param problem Problem for which extreme class cardinalities will be calculated.
#' @return \emph{p} x \emph{2} matrix, where \emph{p} is the number of classes.
#' Value at \code{[h, 1]} is a minimal possible cardinality of class \code{C_h},
#' and value at \code{[h, 2]} is a maximal possible cardinality of class \code{C_h}.
#' @seealso
#' \code{\link{addMinimalClassCardinalities}}
#' \code{\link{addMaximalClassCardinalities}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' extremeClassCardinalities <- calculateExtremeClassCardinalities(problem)
#' @export
calculateExtremeClassCardinalities <- function(problem) {
  if (!checkConsistency(problem))
    stop("Model infeasible.")
  
  model <- buildModel(problem, FALSE)
  
  if (ncol(model$constraints$lhs) == model$firstThresholdIndex + problem$nrClasses - 2) {
    model <- extendModelWithAssignmentVariables(model)
  }

  result <- matrix(nrow = problem$nrClasses, ncol = 2)
  
  for(class in 1:problem$nrClasses) {
    result[class, 1] <- calculateExtremeClassCardinality(model, class, FALSE)
    if (RORUTADIS_VERBOSE) print (paste("minimal cardinality of class ", class, "equals", rel[class, 1]))
    
    result[class, 2] <- calculateExtremeClassCardinality(model, class, TRUE)
    if (RORUTADIS_VERBOSE) print (paste("maximal cardinality of class ", class, "equals", rel[class, 2]))
  }
  
  return (result)
}

#' Merge different assignments
#'
#' This function allows to merge different assignments, e.g. from various
#' decision makers (group result, group assignment).
#' There are four types of group assignments:
#' \itemize{
#' \item \strong{P}ossible \strong{P}ossible -
#' alternative \emph{a_i} is \strong{possibly} in class \emph{C_h}
#' \strong{for at least one} decision maker,
#' \item \strong{P}ossible \strong{N}ecessary - 
#' alternative \emph{a_i} is \strong{possibly} in class \emph{C_h}
#' \strong{for all} decision makers,
#' \item \strong{N}ecessary \strong{P}ossible - 
#' alternative \emph{a_i} is \strong{necessarily} in class \emph{C_h}
#' \strong{for at least one} decision maker,
#' \item \strong{N}ecessary \strong{N}ecessary - 
#' alternative \emph{a_i} is \strong{necessarily} in class \emph{C_h}
#' \strong{for all} decision makers.
#' }
#' The first possible-necessary parameter depends on decision makers
#' assignments computed earlier, and the second is define as function parameter.
#' 
#' @param assignmentList List of assignment matrices (results of calling
#' \code{\link{calculateAssignments}} function).
#' @param necessary Whether necessary or possible merging.
#' @return \emph{n} x \emph{p} logical matrix, where each row represents one
#' of \emph{n} alternatives and each column represents one of \emph{p} classes.
#' Element \code{[i, h]} is \code{TRUE} if alternative \code{a_i} can be assigned
#' to class \code{C_h}.
#' @seealso
#' \code{\link{calculateAssignments}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' DM1Problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' DM2Problem <- addAssignmentsLB(problem, c(2, 2), c(4, 2))
#' 
#' necessary <- FALSE
#' assignmentList <- list()
#' assignmentList[[1]] <- calculateAssignments(DM1Problem, necessary)
#' assignmentList[[2]] <- calculateAssignments(DM2Problem, necessary)
#' 
#' # generate possible - necessary assignments
#' PNAssignments <- mergeAssignments(assignmentList, TRUE)
#' @export
mergeAssignments <- function(assignmentList, necessary) {
  stopifnot(is.list(assignmentList))
  stopifnot(length(assignmentList) > 0)
  stopifnot(is.logical(necessary))
  
  if (length(assignmentList) == 1)
    return (assignmentList[1])
  
  result <- assignmentList[[1]]
  
  if (necessary) {
    for (k in 2:length(assignmentList)) {
      for (h in 1:ncol(result)) {
        for (i in 1:nrow(result)) {
          result[i, h] <- result[i, h] && assignmentList[[k]][i, h]
        }
      }
    }
  }
  else {
    for (k in 2:length(assignmentList)) {
      for (h in 1:ncol(result)) {
        for (i in 1:nrow(result)) {
          result[i, h] <- result[i, h] || assignmentList[[k]][i, h]
        }
      }
    }
  }
  
  return (result)
}

#' Find representative utility function
#'
#' This function finds a representative utility function for a problem.
#' 
#' @param problem Problem to investigate.
#' @param mode An integer that represents a method of a computing representative
#' utility function:
#' \itemize{
#' \item \code{0} - iterative mode,
#' \item \code{1} - compromise mode.
#' }
#' @param relation A matrix of assignment pairwise comparisons (see
#' \code{\link{compareAssignments}}). If the parameter is \code{NULL}, the relation
#' will be computed.
#' @return List with named elements:
#' \itemize{
#' \item \code{vf} - list of 2-column matrices with marginal value functions (characteristic point in rows),
#' \item \code{thresholds},
#' \item \code{assignments},
#' \item \code{alternativeValues},
#' \item \code{epsilon}.
#' }
#' \code{NULL} is returned if representative function cannot be found.
#' @seealso
#' \code{\link{plotVF}}
#' \code{\link{plotComprehensiveValue}}
#' \code{\link{findSimpleFunction}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' representativeFunction <- findRepresentativeFunction(problem, 0)
#' assignments <- representativeFunction$assignments
#' @export
findRepresentativeFunction <- function(problem, mode, relation = NULL) {
  stopifnot(mode == 0 || mode == 1)
  
  solution <- NULL
  
  if (is.null(relation)) {
    relation <- compareAssignments(problem)
  }
  
  model <- buildModel(problem, F)
  nrAlternatives <- nrow(problem$perf)
  
  if (mode == 0) { # iterative mode    
    deltaIndex <- ncol(model$constraints$lhs) + 1
    gammaIndex <- ncol(model$constraints$lhs) + 3
    model$constraints <- addVarialbesToModel(model$constraints, c("C", "C", "C", "C"))
    nrVariables <- ncol(model$constraints$lhs)
    
    for (i in 1:(nrAlternatives - 1)) {
      for (j in (i + 1):nrAlternatives) {
        if (relation[i, j] && !relation[j, i]) {
          #U(a_i) - U(a_j) >= delta
          valueDifferenceLhs <- ua(i, nrVariables, model$perfToModelVariables) -
            ua(j, nrVariables, model$perfToModelVariables)
          valueDifferenceLhs[deltaIndex] <- -1
          valueDifferenceLhs[deltaIndex + 1] <- 1
          model$constraints <- combineConstraints(model$constraints,
                                                  list(lhs = valueDifferenceLhs, dir = ">=", rhs = 0))
        } else if (!relation[i, j] && relation[j, i]) {
          #U(a_j) - U(a_i) >= delta
          valueDifferenceLhs <- ua(j, nrVariables, model$perfToModelVariables) -
            ua(i, nrVariables, model$perfToModelVariables)
          valueDifferenceLhs[deltaIndex] <- -1
          valueDifferenceLhs[deltaIndex + 1] <- 1
          model$constraints <- combineConstraints(model$constraints,
                                                  list(lhs = valueDifferenceLhs, dir = ">=", rhs = 0))
        }
      }
    }
    
    obj <- rep(0, nrVariables)
    obj[deltaIndex] <- 1  
    obj[deltaIndex + 1] <- -1 
    solution <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                               max = TRUE, types = model$constraints$types)
    
    if (solution$status != 0) {
      return (NULL)
    }
    
    if (solution$optimum < 0) {
      warning(paste("Solution was found, but the optimization target is negative (", solution$optimum, ").", sep = ""))
    }
    
    fixedDeltaConstr <- rep(0, nrVariables)
    fixedDeltaConstr[deltaIndex] <- 1
    fixedDeltaConstr[deltaIndex + 1] <- -1
    model$constraints <- combineConstraints(model$constraints, list(lhs = fixedDeltaConstr,
                                                                    dir = "==",
                                                                    rhs = solution$optimum))
    
    for (i in 1:(nrAlternatives - 1)) {
      for (j in (i + 1):nrAlternatives) {
        if (relation[i, j] == relation[j, i]) {
          #U(a_i) - U(a_j) <= gamma
          #U(a_j) - U(a_i) <= gamma
          ijValueDifferenceLhs <- ua(i, nrVariables, model$perfToModelVariables) -
            ua(j, nrVariables, model$perfToModelVariables)
          ijValueDifferenceLhs[gammaIndex] <- -1
          ijValueDifferenceLhs[gammaIndex + 1] <- 1
          
          jiValueDifferenceLhs <- ua(j, nrVariables, model$perfToModelVariables) -
            ua(i, nrVariables, model$perfToModelVariables)
          jiValueDifferenceLhs[gammaIndex] <- -1
          jiValueDifferenceLhs[gammaIndex + 1] <- 1
          
          model$constraints <- combineConstraints(model$constraints,
                                                  list(lhs = ijValueDifferenceLhs,
                                                       dir = "<=",
                                                       rhs = 0),
                                                  list(lhs = jiValueDifferenceLhs,
                                                       dir = "<=",
                                                       rhs = 0))
        }
      }
    }
    
    obj <- rep(0, nrVariables)
    obj[gammaIndex] <- 1  
    obj[gammaIndex + 1] <- -1 
    solution <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                               max = FALSE, types = model$constraints$types)
    
    if (solution$status != 0) {
      return (NULL)
    }
    
    if (solution$optimum < 0) {
      warning(paste("Solution was found, but the optimization target is negative (", solution$optimum, ").", sep = ""))
    }
  } else if (mode == 1) { # compromise mode    
    deltaIndex <- ncol(model$constraints$lhs) + 1
    model$constraints <- addVarialbesToModel(model$constraints, c("C", "C"))
    nrVariables <- ncol(model$constraints$lhs)
    different <- c()
    similar <- c()
    
    for (i in 1:(nrAlternatives - 1)) {
      for (j in (i + 1):nrAlternatives) {
        if (relation[i, j] && !relation[j, i]) {
          different <- c(different, i, j)
        } else if (!relation[i, j] && relation[j, i]) {
          different <- c(different, j, i)
        } else {
          similar <- c(similar, i, j)
        }
      }
    }
    
    if (length(similar) == 0 || length(different) == 0) {
      stop("Compromise mode (1) does not work in this case. Try with another mode.")
    }
    
    for (i in seq(1, length(different), by = 2)) {
      for (j in seq(1, length(similar), by = 2)) {
        valueDifferenceLhs <- ua(different[i], nrVariables, model$perfToModelVariables) -
          ua(different[i + 1], nrVariables, model$perfToModelVariables) +
          ua(similar[j], nrVariables, model$perfToModelVariables) -
          ua(similar[j + 1], nrVariables, model$perfToModelVariables)
        valueDifferenceLhs[deltaIndex] <- -1
        valueDifferenceLhs[deltaIndex + 1] <- 1
        model$constraints <- combineConstraints(model$constraints,
                                                list(lhs = valueDifferenceLhs, dir = ">=", rhs = 0))
      }
    }
    
    obj <- rep(0, nrVariables)
    obj[deltaIndex] <- 1  
    obj[deltaIndex + 1] <- -1 
    solution <- Rglpk_solve_LP(obj, model$constraints$lhs, model$constraints$dir, model$constraints$rhs,
                               max = TRUE, types = model$constraints$types)
    
    if (solution$status != 0) {
      return (NULL)
    }
    
    if (solution$optimum < 0) {
      warning(paste("Solution was found, but the optimization target is negative (", solution$optimum, ").", sep = ""))
    }
  }
  
  return (toSolution(model, solution$solution))
}

#' Find one value function
#'
#' This function finds single value function that is consistent with provided preferece information.
#' Search is done by epsilon maximization.
#' 
#' @param problem Problem
#' @return List with named elements:
#' \itemize{
#' \item \code{vf} - list of 2-column matrices with marginal value functions (characteristic point in rows),
#' \item \code{thresholds},
#' \item \code{assignments},
#' \item \code{alternativeValues},
#' \item \code{epsilon}.
#' }
#' @seealso
#' \code{\link{plotVF}}
#' \code{\link{plotComprehensiveValue}}
#' \code{\link{findRepresentativeFunction}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' simpleFunction <- findSimpleFunction(problem)
#' @export
findSimpleFunction <- function(problem) {
  model <- buildModel(problem, TRUE)
  solution <- maximizeEpsilon(model)
  
  if (solution$status == 0 && solution$optimum >= RORUTADIS_MINEPS) {
    return (toSolution(model, solution$solution))
  }
  
  return (NULL)
}

