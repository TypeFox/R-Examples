#### HELPERS

isSuperset <- function(set, subset) {
  for (i in subset)
    if (!(i %in% set))
      return (FALSE)
  return (TRUE)
}

#### EXPLANATIONS

#' Explain assignment
#'
#' This function allows to obtain explanation of an alternative assignment to
#' a specific class interval or one class in case if assignment is necessary.
#' The function returns all preferential reducts for an assignment relation.
#'
#' @param alternative Index of an alternative.
#' @param classInterval Two-element vector \code{c(l, u)} that represents
#' an assignment of \code{alternative} to class interval \code{[C_l, C_u]}
#' (\code{l <= u}). 
#' @param problem Problem for which computations will be performed.
#' @return List of all preferential reducts for an assignment relation.
#' If the assignment is not influenced by restrictions then empty list will be returned.
#' Each element of the list is a preferential reduct represented as a vector
#' of restriction indices. To identify preferential core use
#' \code{\link{getPreferentialCore}}.
#' To find out about restrictions by their indices use \code{\link{getRestrictions}}.
#' @seealso
#' \code{\link{getPreferentialCore}}
#' \code{\link{getRestrictions}}
#' \code{\link{calculateAssignments}}
#' @importFrom utils combn
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.5), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' possibleAssignments <- calculateAssignments(problem, FALSE)
#' alternative <- 4
#' assignment <- c(min(which(possibleAssignments[alternative, ])),
#'                max(which(possibleAssignments[alternative, ])))
#'                
#' preferentialReducts <- explainAssignment(alternative,
#'    assignment, problem)
#' preferentialCore <- getPreferentialCore(preferentialReducts)
#' coreRestrictions <- getRestrictions(problem, preferentialCore)
#' @export
explainAssignment <- function(alternative, classInterval, problem) {
  stopifnot(is.vector(classInterval))
  stopifnot(length(classInterval) == 2)
  stopifnot(classInterval[1] <= classInterval[2])
  
  # it is not good to have classInterval as parameter
  # but it was left for compatibility with older versions for now
  # however, it is checked for corectness:
  
  model <- buildModel(problem, T)
  assignments <- sapply(seq_len(problem$nrClasses), function(h) { checkRelation(model, alternative, h, F) })
  assignmentInterval <- c(min(which(assignments)), max(which(assignments)))
  stopifnot(assignmentInterval[1] == classInterval[1])
  stopifnot(assignmentInterval[2] == classInterval[2])
  stopifnot(all(assignments[assignmentInterval[1]:assignmentInterval[2]]))
  
  fromClass <- assignmentInterval[1]
  toClass <- assignmentInterval[2]
  
  stopifnot(fromClass >= 1)
  stopifnot(toClass >= 1)
  stopifnot(fromClass <= problem$nrClasses)
  stopifnot(toClass <= problem$nrClasses)
  
  if (!isModelConsistent(model)) {
    stop("The problem is inconsistent.")
  }
  
  if (fromClass == 1 && toClass < problem$nrClasses) {
    model$constraints <- combineConstraints(model$constraints,
                                            buildLBAssignmentsConstraint(alternative, toClass + 1, model))
  } else if (fromClass > 1 && toClass == problem$nrClasses) {
    model$constraints <- combineConstraints(model$constraints,
                                            buildUBAssignmentsConstraint(alternative, fromClass - 1, model))
  } else if (fromClass > 1 && toClass < problem$nrClasses) {
    model$constraints <- addVarialbesToModel(model$constraints, c("B", "B"))
    nrVariables <- ncol(model$constraints$lhs)
    
    constr1 <- buildLBAssignmentsConstraint(alternative, toClass + 1, model)
    constr2 <- buildUBAssignmentsConstraint(alternative, fromClass - 1, model)
    constrRel <- list(lhs = rep(0, nrVariables), dir = "==", rhs = 1)
    
    constr1$lhs[nrVariables - 1] <- RORUTADIS_BIGM
    constr2$lhs[nrVariables] <- -RORUTADIS_BIGM
    constrRel$lhs[nrVariables - 1] <- 1
    constrRel$lhs[nrVariables] <- 1
    
    model$constraints <- combineConstraints(model$constraints, constr1, constr2, constrRel)
  }
  
  nrPreferenceInformation <- length(model$prefInfoToConstraints)
  preferentialReducts <- list()
  
  if (nrPreferenceInformation > 0) {
    subsets <- lapply(seq_len(nrPreferenceInformation), function(x) combn(nrPreferenceInformation, x))
    
    newModel <- model
    newModel$constraints <- removeConstraints(newModel$constraints, unlist(model$prefInfoToConstraints))
    
    if (!isModelConsistent(newModel)) {
      return (list())
    }
    
    for (i in 1:length(subsets)) {
      for (j in 1:ncol(subsets[[i]])) {
        subset <- subsets[[i]][, j]
        if (subset[1] > 0) {
          newModel <- model
          
          if (length(subset) != nrPreferenceInformation) {
            newModel$constraints <- removeConstraints(newModel$constraints, unlist(model$prefInfoToConstraints[-subset]))
          }
          
          if (!isModelConsistent(newModel)) {
            preferentialReducts[[length(preferentialReducts) + 1]] <- subset
            
            for (k in 1:length(subsets)) {
              for (l in 1:ncol(subsets[[k]])) {
                if (isSuperset(subsets[[k]][, l], subset)) {
                  subsets[[k]][1, l] <- 0
                }
              }
            }
          }
        }
      }
    }
  }
  
  return (preferentialReducts)
}

#' Identify preferential core
#'
#' This function identifies preferential core.
#'
#' @param preferentialReducts List of all preferential reducts 
#' (a result of \code{\link{explainAssignment}}).
#' @return Preferential core as a vector of restriction indices.
#' To find out about restrictions by their indices use \code{\link{getRestrictions}}.
#' @seealso
#' \code{\link{explainAssignment}}
#' \code{\link{getRestrictions}}
#' @examples
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.5), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' possibleAssignments <- calculateAssignments(problem, FALSE)
#' alternative <- 4
#' assignment <- c(min(which(possibleAssignments[alternative, ])),
#'                max(which(possibleAssignments[alternative, ])))
#'                
#' preferentialReducts <- explainAssignment(alternative,
#'    assignment, problem)
#' preferentialCore <- getPreferentialCore(preferentialReducts)
#' coreRestrictions <- getRestrictions(problem, preferentialCore)
#' @export
getPreferentialCore <- function(preferentialReducts) {
  if (length(preferentialReducts) == 0) {
    return (c())
  }
  else if (length(preferentialReducts) == 1) {
    return (preferentialReducts[[1]])
  }
  
  res <- preferentialReducts[[1]]
  
  for (i in 2:length(preferentialReducts)) {
    oldRes <- res
    res <- c()
    for (j in 1:length(oldRes)) {
      if (oldRes[j] %in% preferentialReducts[[i]]) {
        res <- c(res, oldRes[j])
      }
    }
  }
  
  return (res)
}
