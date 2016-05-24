#### INCREMENTAL PROBLEM BUILDING

#' Build a representation of a problem
#'
#' This function creates representation of a given problem for usage
#' in farther computations.
#'
#' @param perf A \emph{n} x \emph{m} performance matrix of \emph{n} alternatives evaluated
#' on \emph{m} criteria.
#' @param nrClasses Number of classes.
#' @param strictVF \code{TRUE} for strictly monotonic marginal value functions,
#' \code{FALSE} for weakly monotonic.
#' @param criteria A vector containing type of each criterion (\code{'g'} - gain, \code{'c'} - cost).
#' @param characteristicPoints A vector of integers that for each criterion contains number of characteristic points
#' or \emph{0} for general marginal value function.
#' @return Representation of a problem as a list with named members.
#' @seealso
#' \code{\link{addAssignmentsLB}}
#' \code{\link{removeAssignmentsLB}}
#' \code{\link{addAssignmentsUB}}
#' \code{\link{removeAssignmentsUB}}
#' \code{\link{addAssignmentPairwiseAtLeastComparisons}}
#' \code{\link{removeAssignmentPairwiseAtLeastComparisons}}
#' \code{\link{addAssignmentPairwiseAtMostComparisons}}
#' \code{\link{removeAssignmentPairwiseAtMostComparisons}}
#' \code{\link{addMinimalClassCardinalities}}
#' \code{\link{removeMinimalClassCardinalities}}
#' \code{\link{addMaximalClassCardinalities}}
#' \code{\link{removeMaximalClassCardinalities}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' @export
buildProblem <- function(perf, nrClasses, strictVF, criteria, characteristicPoints) {
  stopifnot(is.matrix(perf))
  stopifnot(is.vector(characteristicPoints))
  stopifnot(ncol(perf) == length(criteria))
  stopifnot(length(which(criteria != 'g' & criteria != 'c')) == 0)
  stopifnot(ncol(perf) == length(characteristicPoints))
  stopifnot(is.logical(strictVF))
  stopifnot(nrClasses >= 2)
  #stopifnot(all(characteristicPoints >= 0) && all(characteristicPoints != 1))
  
  return (list(perf = perf,
               nrClasses = nrClasses,
               strictVF = strictVF,
               criteria = criteria,
               characteristicPoints = characteristicPoints,
               assignmentsLB = NULL,
               assignmentsUB = NULL,
               assignmentPairwiseAtLeastComparisons = NULL,
               assignmentPairwiseAtMostComparisons = NULL,
               minimalClassCardinalities = NULL,
               maximalClassCardinalities = NULL))
}

###### assignmentsLB

#' Add lower bound of alternative possible assignments
#'
#' This function adds lower bounds of possible assignments to
#' a problem.
#'
#' @param problem Problem to which preference information will be added.
#' @param ... Assignments as two-element vectors.
#' Each vector \code{c(i, j)} represents assignment of an alternative \emph{a_i} 
#' to class at least as good as class \emph{C_j}.
#' @return Problem with added assignment examples.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeAssignmentsLB}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add assignment examples: alternative 1 to class at least as good as class 2
#' # and alternative 2 to class at least as good as class 3
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' @export
addAssignmentsLB <- function(problem, ...) {
  assignments <- list(...)
  
  for (assignment in assignments) {
    stopifnot(length(assignment) == 2)
    stopifnot(assignment[1] > 0)
    stopifnot(assignment[2] > 0)
    stopifnot(assignment[1] <= nrow(problem$perf))
    stopifnot(assignment[2] <= problem$nrClasses)
    
    if (is.null(problem$assignmentsLB)) {
      problem$assignmentsLB <- matrix(assignment, ncol = 2)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$assignmentsLB)) {
        if (problem$assignmentsLB[i, 1] == assignment[1]) {
          problem$assignmentsLB[i, 2] <- assignment[2]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$assignmentsLB <- rbind(problem$assignmentsLB,
                                       assignment, deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove lower bound of alternative possible assignments
#'
#' This function removes lower bounds of possible assignments from
#' a problem.
#'
#' @param problem Problem from which preference information will be removed.
#' @param ... Assignments as two-element vectors and/or integers.
#' Each argument represents assignment to remove. If  \code{c(i, j)} vector was
#' provided an assignment of an alternative \emph{a_i} 
#' to at least class \emph{C_j} will be removed. In case where single value \code{i} was
#' given an assignment of an alternative \emph{a_i} will be removed regardless of class.
#' If a specific assignment was not found nothing will happen.
#' 
#' @return Problem with removed assignment examples.
#' 
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add assignment examples: alternative 1 at least to class 2
#' # alternative 2 at least to class 3
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' # and remove the assignments
#' problem <- removeAssignmentsLB(problem, c(1, 2), 2)
#' @export
removeAssignmentsLB <- function(problem, ...) {
  assignments <- list(...)
  
  for (assignment in assignments) {
    stopifnot(length(assignment) == 1 || length(assignment) == 2)
    
    if (!is.null(problem$assignmentsLB)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$assignmentsLB)) {
        if (length(assignment) == 1 &&
              problem$assignmentsLB[i, 1] != assignment[1]) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentsLB[i, ], deparse.level = 0)
        }
        else if (length(assignment) == 2 &&
                   (problem$assignmentsLB[i, 1] != assignment[1] ||
                      problem$assignmentsLB[i, 2] != assignment[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentsLB[i, ], deparse.level = 0)
        }
      }
      
      problem$assignmentsLB <- tmpRestrictions
    }
  }
  
  return (problem)
}

###### assignmentsUB

#' Add upper bound of alternative possible assignments
#'
#' This function adds upper bounds of possible assignments to a problem.
#'
#' @param problem Problem to which preference information will be added.
#' @param ... Assignments as two-element vectors.
#' Each vector \code{c(i, j)} represents assignment of an alternative \emph{a_i} 
#' to at most class as good as \emph{C_j}.
#' @return Problem with added assignment examples.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeAssignmentsUB}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add assignment examples: alternative 3 at most to class as good as class 1
#' # and alternative 4 to class at most as good as class 2
#' problem <- addAssignmentsUB(problem, c(3, 1), c(4, 2))
#' @export
addAssignmentsUB <- function(problem, ...) {
  assignments <- list(...)
  
  for (assignment in assignments) {
    stopifnot(length(assignment) == 2)
    stopifnot(assignment[1] > 0)
    stopifnot(assignment[2] > 0)
    stopifnot(assignment[1] <= nrow(problem$perf))
    stopifnot(assignment[2] <= problem$nrClasses)
    
    if (is.null(problem$assignmentsUB)) {
      problem$assignmentsUB <- matrix(assignment, ncol = 2)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$assignmentsUB)) {
        if (problem$assignmentsUB[i, 1] == assignment[1]) {
          problem$assignmentsUB[i, 2] <- assignment[2]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$assignmentsUB <- rbind(problem$assignmentsUB,
                                       assignment, deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove upper bound of alternative possible assignments
#'
#' This function removes upper bounds of possible assignments from
#' a problem.
#'
#' @param problem Problem from which preference information will be removed.
#' @param ... Assignments as two-element vectors and/or integers.
#' Each argument represents assignment to remove. If  \code{c(i, j)} vector was
#' provided an assignment of an alternative \emph{a_i} 
#' to at most class \emph{C_j} will be removed. In case where single value \code{i} was
#' given an assignment of an alternative \emph{a_i} will be removed regardless of class.
#' If a specific assignment was not found nothing will happen.
#' 
#' @return Problem with removed assignment examples.
#' 
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add assignment examples: alternative 1 at least to class 2
#' # alternative 2 at least to class 3
#' problem <- addAssignmentsLB(problem, c(1, 2), c(2, 3))
#' 
#' # and remove the assignments
#' problem <- removeAssignmentsLB(problem, c(1, 2), 2)
#' @export
removeAssignmentsUB <- function(problem, ...) {
  assignments <- list(...)
  
  for (assignment in assignments) {
    stopifnot(length(assignment) == 1 || length(assignment) == 2)
    
    if (!is.null(problem$assignmentsUB)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$assignmentsUB)) {
        if (length(assignment) == 1 &&
              problem$assignmentsUB[i, 1] != assignment[1]) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentsUB[i, ], deparse.level = 0)
        }
        else if (length(assignment) == 2 &&
                   (problem$assignmentsUB[i, 1] != assignment[1] ||
                      problem$assignmentsUB[i, 2] != assignment[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentsUB[i, ], deparse.level = 0)
        }
      }
      
      problem$assignmentsUB <- tmpRestrictions
    }
  }
  
  return (problem)
}

###### assignmentPairwiseAtLeastComparisons

#' Add assignment pairwise \emph{at least} comparisons
#'
#' The comparison of a pair of alternatives may indicate that \emph{a_i} should
#' be assigned to a class at least as good as class of \emph{a_j} or at least
#' better by \emph{k} classes. The function \code{assignmentPairwiseAtLeastComparisons}
#' allows to define such pairwise comparisons.
#' 
#' @param problem Problem to which preference information will be added.
#' @param ... Comparisons as three-element vectors.
#' Each vector \code{c(i, j, k)} represents a single assignment comparison:
#' alternative \emph{a_i} has to be assigned to class at least better by
#' \emph{k} classes then class of \emph{a_j}.
#' @return Problem with added comparisons.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeAssignmentPairwiseAtLeastComparisons}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add comparisons:
#' # alternative 2 to class at least as good as class of alternative 1
#' # alternative 4 to class at least better by 1 class then class
#' # of alternative 3
#' problem <- addAssignmentPairwiseAtLeastComparisons(problem, 
#'    c(4, 3, 1), c(2, 1, 0))
#' @export
addAssignmentPairwiseAtLeastComparisons <- function(problem, ...) {
  comparisons <- list(...)
  
  for (comparison in comparisons) {
    stopifnot(length(comparison) == 3)
    stopifnot(comparison[1] > 0)
    stopifnot(comparison[2] > 0)
    stopifnot(comparison[3] >= 0)
    stopifnot(comparison[1] <= nrow(problem$perf))
    stopifnot(comparison[2] <= nrow(problem$perf))
    stopifnot(comparison[3] <= problem$nrClasses - 1)
    
    if (is.null(problem$assignmentPairwiseAtLeastComparisons)) {
      problem$assignmentPairwiseAtLeastComparisons <- matrix(comparison, ncol = 3)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$assignmentPairwiseAtLeastComparisons)) {
        if (problem$assignmentPairwiseAtLeastComparisons[i, 1] == comparison[1] &&
              problem$assignmentPairwiseAtLeastComparisons[i, 2] == comparison[2]) {
          problem$assignmentPairwiseAtLeastComparisons[i, 3] <- comparison[3]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$assignmentPairwiseAtLeastComparisons <-rbind(problem$assignmentPairwiseAtLeastComparisons,
                                                             comparison,
                                                             deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove assignment pairwise \emph{at least} comparisons
#'
#' This function removes pairwise \emph{at least} comparisons. For more
#' information see \code{addPairwiseAtLeastComparisons}.
#' 
#' @param problem Problem from which preference information will be removed
#' @param ... Comparisons as three-element vectors and/or two-element vectors.
#' Each argument represents comparison to remove. If \code{c(i, j, k)} vector was
#' provided a corresponding comparison will be removed. In case where two-element
#' vector \code{c(i,j)} was given a comparison of an alternative \emph{a_i} with
#' \emph{a_j} will be removed regardless of value of \emph{k}.
#' If a specific comparison was not found nothing will happen.
#' 
#' @return Problem with removed comparisons.
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add comparisons:
#' # alternative 2 to class at least as good as class of alternative 1
#' # alternative 4 to class at least better by 1 class then class
#' # of alternative 3
#' problem <- addAssignmentPairwiseAtLeastComparisons(problem, 
#'    c(4, 3, 1), c(2, 1, 0))
#' # remove comparison between alternative 4 and 3
#' problem <- removeAssignmentPairwiseAtLeastComparisons(problem, c(4, 3))
#' @export
removeAssignmentPairwiseAtLeastComparisons <- function(problem, ...) {
  comparisons <- list(...)
  
  for (comparison in comparisons) {
    stopifnot(length(comparison) == 2 || length(comparison) == 3)
    
    if (!is.null(problem$assignmentPairwiseAtLeastComparisons)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$assignmentPairwiseAtLeastComparisons)) {
        if (length(comparison) == 2 &&
              (problem$assignmentPairwiseAtLeastComparisons[i, 1] != comparison[1] ||
                 problem$assignmentPairwiseAtLeastComparisons[i, 2] != comparison[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentPairwiseAtLeastComparisons[i, ],
                                   deparse.level = 0)
        }
        else if (length(comparison) == 3 &&
                   (problem$assignmentPairwiseAtLeastComparisons[i, 1] != comparison[1] ||
                      problem$assignmentPairwiseAtLeastComparisons[i, 2] != comparison[2] ||
                      problem$assignmentPairwiseAtLeastComparisons[i, 3] != comparison[3])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentPairwiseAtLeastComparisons[i, ],
                                   deparse.level = 0)
        }
      }
      
      problem$assignmentPairwiseAtLeastComparisons <- tmpRestrictions
    }
  }
  
  return (problem)
}

###### assignmentPairwiseAtMostComparisons

#' Add assignment pairwise \emph{at most} comparisons
#'
#' The comparison of a pair of alternatives may indicate that alternative
#' \emph{a_i} should be assigned to a class at most better by \emph{k} classes
#' then class of \emph{a_j}. The function \code{assignmentPairwiseAtMostComparisons}
#' allows to define such pairwise comparisons.
#' 
#' @param problem Problem to which preference information will be added.
#' @param ... Comparisons as three-element vectors.
#' Each vector \code{c(i, j, k)} represents a single assignment comparison:
#' alternative \emph{a_i} has to be assigned to class at most better by
#' \emph{k} classes then class of \emph{a_j}.
#' @return Problem with added comparisons.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeAssignmentPairwiseAtMostComparisons}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add comparison:
#' # alternative 4 to class at most better by 1 class then class
#' # of alternative 3
#' problem <- addAssignmentPairwiseAtMostComparisons(problem, c(4, 3, 1))
#' @export
addAssignmentPairwiseAtMostComparisons <- function(problem, ...) {
  comparisons <- list(...)
  
  for (comparison in comparisons) {
    stopifnot(length(comparison) == 3)
    stopifnot(comparison[1] > 0)
    stopifnot(comparison[2] > 0)
    stopifnot(comparison[3] >= 0)
    stopifnot(comparison[1] <= nrow(problem$perf))
    stopifnot(comparison[2] <= nrow(problem$perf))
    stopifnot(comparison[3] <= problem$nrClasses - 1)
    
    if (is.null(problem$assignmentPairwiseAtMostComparisons)) {
      problem$assignmentPairwiseAtMostComparisons <- matrix(comparison, ncol = 3)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$assignmentPairwiseAtMostComparisons)) {
        if (problem$assignmentPairwiseAtMostComparisons[i, 1] == comparison[1] &&
              problem$assignmentPairwiseAtMostComparisons[i, 2] == comparison[2]) {
          problem$assignmentPairwiseAtMostComparisons[i, 3] <- comparison[3]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$assignmentPairwiseAtMostComparisons <- rbind(problem$assignmentPairwiseAtMostComparisons,
                                                             comparison,
                                                             deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove assignment pairwise \emph{at most} comparisons
#'
#' This function removes pairwise \emph{at most} comparisons. For more
#' information see \code{addPairwiseAtMostComparisons}.
#' 
#' @param problem Problem from which preference information will be removed
#' @param ... Comparisons as three-element vectors and/or two-element vectors.
#' Each argument represents comparison to remove. If \code{c(i, j, k)} vector was
#' provided a corresponding comparison will be removed. In case where two-element
#' vector \code{c(i,j)} was given a comparison of an alternative \emph{a_i} with
#' \emph{a_j} will be removed regardless of value of \emph{k}.
#' If a specific comparison was not found nothing will happen.
#' 
#' @return Problem with removed comparisons.
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # add comparison:
#' # alternative 4 to class at most better by 1 class then class
#' # of alternative 3
#' problem <- addAssignmentPairwiseAtMostComparisons(problem, c(4, 3, 1))
#' # remove comparison between alternative 4 and 3
#' problem <- removeAssignmentPairwiseAtMostComparisons(problem, c(4, 3))
#' @export
removeAssignmentPairwiseAtMostComparisons <- function(problem, ...) {
  comparisons <- list(...)
  
  for (comparison in comparisons) {
    stopifnot(length(comparison) == 2 || length(comparison) == 3)
    
    if (!is.null(problem$assignmentPairwiseAtMostComparisons)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$assignmentPairwiseAtMostComparisons)) {
        if (length(comparison) == 2 &&
              (problem$assignmentPairwiseAtMostComparisons[i, 1] != comparison[1] ||
                 problem$assignmentPairwiseAtMostComparisons[i, 2] != comparison[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentPairwiseAtMostComparisons[i, ],
                                   deparse.level = 0)
        }
        else if (length(comparison) == 3 &&
                   (problem$assignmentPairwiseAtMostComparisons[i, 1] != comparison[1] ||
                      problem$assignmentPairwiseAtMostComparisons[i, 2] != comparison[2] ||
                      problem$assignmentPairwiseAtMostComparisons[i, 3] != comparison[3])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$assignmentPairwiseAtMostComparisons[i, ],
                                   deparse.level = 0)
        }
      }
      
      problem$assignmentPairwiseAtMostComparisons <- tmpRestrictions
    }
  }
  
  return (problem)
}

###### minimalClassCardinalities

#' Add minimal class cardinality restrictions
#'
#' This function allows to define minimal cardinality of particular classes.
#'
#' @param problem Problem to which preference information will be added.
#' @param ... Minimal cardinalities as two-element vectors \code{c(i, j)}, where
#' \emph{j} is a minimal cardinality of class \emph{C_i}.
#' @return Problem with added preference information.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeMinimalClassCardinalities}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # set minimal class cardinalities:
#' # at least one alternative has to be assigned to class 2
#' # and at least one alternative has to be assigned to class 3
#' problem <- addMinimalClassCardinalities(problem, c(2, 1), c(3, 1))
#' @export
addMinimalClassCardinalities <- function(problem, ...) {
  restrictions <- list(...)
  
  for (cardinalityRestriction in restrictions) {
    stopifnot(length(cardinalityRestriction) == 2)
    stopifnot(cardinalityRestriction[1] > 0)
    stopifnot(cardinalityRestriction[2] > 0)
    stopifnot(cardinalityRestriction[1] <= problem$nrClasses)
    
    if (is.null(problem$minimalClassCardinalities)) {
      problem$minimalClassCardinalities <- matrix(cardinalityRestriction, ncol = 2)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$minimalClassCardinalities)) {
        if (problem$minimalClassCardinalities[i, 1] == cardinalityRestriction[1]) {
          problem$minimalClassCardinalities[i, 2] <- cardinalityRestriction[2]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$minimalClassCardinalities <- rbind(problem$minimalClassCardinalities,
                                                     cardinalityRestriction,
                                                     deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove minimal class cardinality restrictions
#'
#' This function allows to remove defined minimal cardinality of particular classes.
#'
#' @param problem Problem from which preference information will be removed.
#' @param ... Two-element vectors and/or integers.
#' Each argument represents restriction to remove. If  \code{c(i, j)} vector was
#' provided then defined minimal cardinality \emph{j} for class \emph{C_i} will
#' be removed. In case where single value \code{i} was given a restriction for
#' class \emph{a_i} will be removed regardless of minimal cardinality value.
#' If a specific restriction was not found nothing will happen.
#' 
#' @return Problem with removed preference information.
#' 
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # set minimal class cardinalities:
#' # at least one alternative has to be assigned to class 2
#' # and at least one alternative has to be assigned to class 3
#' problem <- addMinimalClassCardinalities(problem, c(2, 1), c(3, 1))
#' # remove defined restriction for class 2
#' problem <- removeMinimalClassCardinalities(problem, 2)
#' @export
removeMinimalClassCardinalities <- function(problem, ...) {
  restrictions <- list(...)
  
  for (cardinalityRestriction in restrictions) {
    stopifnot(length(cardinalityRestriction) == 1 || length(cardinalityRestriction) == 2)
    
    if (!is.null(problem$minimalClassCardinalities)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$minimalClassCardinalities)) {
        if (length(cardinalityRestriction) == 1 &&
              problem$minimalClassCardinalities[i, 1] != cardinalityRestriction[1]) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$minimalClassCardinalities[i, ],
                                   deparse.level = 0)
        }
        else if (length(cardinalityRestriction) == 2 &&
                   (problem$minimalClassCardinalities[i, 1] != cardinalityRestriction[1] ||
                      problem$minimalClassCardinalities[i, 2] != cardinalityRestriction[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$minimalClassCardinalities[i, ],
                                   deparse.level = 0)
        }
      }
      
      problem$minimalClassCardinalities <- tmpRestrictions
    }
  }
  
  return (problem)
}

###### maximalClassCardinalities

#' Add maximal class cardinality restrictions
#'
#' This function allows to define maximal cardinality of particular classes.
#'
#' @param problem Problem to which preference information will be added.
#' @param ... Minimal cardinalities as two-element vectors \code{c(i, j)}, where
#' \emph{j} is a maximal cardinality of class \emph{C_i}.
#' @return Problem with added preference information.
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{removeMaximalClassCardinalities}}
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # set maximal class cardinalities:
#' # at most two alternatives could be assigned to class 2
#' # and at most one alternative could be assigned to class 3
#' problem <- addMaximalClassCardinalities(problem, c(2, 2), c(3, 1))
#' @export
addMaximalClassCardinalities <- function(problem, ...) {
  restrictions <- list(...)
  
  for (cardinalityRestriction in restrictions) {
    stopifnot(length(cardinalityRestriction) == 2)
    stopifnot(cardinalityRestriction[1] > 0)
    stopifnot(cardinalityRestriction[2] > 0)
    stopifnot(cardinalityRestriction[1] <= problem$nrClasses)
    
    if (is.null(problem$maximalClassCardinalities)) {
      problem$maximalClassCardinalities <- matrix(cardinalityRestriction, ncol = 2)
    }
    else {
      found <- FALSE
      for (i in 1:nrow(problem$maximalClassCardinalities)) {
        if (problem$maximalClassCardinalities[i, 1] == cardinalityRestriction[1]) {
          problem$maximalClassCardinalities[i, 2] <- cardinalityRestriction[2]
          found <- TRUE
          break
        }
      }
      if (!found)
        problem$maximalClassCardinalities <- rbind(problem$maximalClassCardinalities,
                                                     cardinalityRestriction,
                                                     deparse.level = 0)
    }
  }
  
  return (problem)
}

#' Remove maximal class cardinality restrictions
#'
#' This function allows to remove defined maximal cardinality of particular classes.
#'
#' @param problem Problem from which preference information will be removed.
#' @param ... Two-element vectors and/or integers.
#' Each argument represents restriction to remove. If  \code{c(i, j)} vector was
#' provided then defined maximal cardinality \emph{j} for class \emph{C_i} will
#' be removed. In case where single value \code{i} was given, a restriction for
#' class \emph{a_i} will be removed regardless of maximal cardinality value.
#' If a specific restriction was not found nothing will happen.
#' 
#' @return Problem with removed preference information.
#' 
#' @examples
#' # 4 alternatives, 2 gain criteria, 3 classes, monotonously increasing
#' # and general marginal value functions
#' perf <- matrix(c(5, 2, 1, 7, 0.5, 0.9, 0.4, 0.4), ncol = 2)
#' problem <- buildProblem(perf, 3, FALSE, c('g', 'g'), c(0, 0))
#' 
#' # set maximal class cardinalities:
#' # at most two alternatives could be assigned to class 2
#' # and at most one alternative could be assigned to class 3
#' problem <- addMaximalClassCardinalities(problem, c(2, 2), c(3, 1))
#' # remove defined restriction for class 2
#' problem <- removeMaximalClassCardinalities(problem, 2)
#' @export
removeMaximalClassCardinalities <- function(problem, ...) {
  restrictions <- list(...)
  
  for (cardinalityRestriction in restrictions) {
    stopifnot(length(cardinalityRestriction) == 1 || length(cardinalityRestriction) == 2)
    
    if (!is.null(problem$maximalClassCardinalities)) {
      tmpRestrictions <- NULL
      
      for (i in 1:nrow(problem$maximalClassCardinalities)) {
        if (length(cardinalityRestriction) == 1 &&
              problem$maximalClassCardinalities[i, 1] != cardinalityRestriction[1]) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$maximalClassCardinalities[i, ],
                                   deparse.level = 0)
        }
        else if (length(cardinalityRestriction) == 2 &&
                   (problem$maximalClassCardinalities[i, 1] != cardinalityRestriction[1] ||
                      problem$maximalClassCardinalities[i, 2] != cardinalityRestriction[2])) {
          tmpRestrictions <- rbind(tmpRestrictions,
                                   problem$maximalClassCardinalities[i, ],
                                   deparse.level = 0)
        }
      }
      
      problem$maximalClassCardinalities <- tmpRestrictions
    }
  }
  
  return (problem)
}

#### GETTING RESTRICTIONS BY INDICES

#' Get restrictions by indices
#'
#' This function gets restrictions by indices. 
#'
#' @param problem Problem whose restrictions will be searched.
#' @param indices A vector of restriction indices (eg. a result of calling
#' \code{\link{getPreferentialCore}}.) Incorrect indices are skipped.
#' @return List with named elements. Each element is a matrix which contains set
#' of restrictions of same type.
#' @seealso
#' \code{\link{getPreferentialCore}}
#' \code{\link{explainAssignment}}
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
getRestrictions <- function(problem, indices) {
  res <- list(assignmentsLB = matrix(ncol = 2, nrow = 0),
              assignmentsUB = matrix(ncol = 2, nrow = 0),
              assignmentPairwiseAtLeastComparisons = matrix(ncol = 3, nrow = 0),
              assignmentPairwiseAtMostComparisons = matrix(ncol = 3, nrow = 0),
              minimalClassCardinalities = matrix(ncol = 2, nrow = 0),
              maximalClassCardinalities = matrix(ncol = 2, nrow = 0))
  index <- 1
  
  if (!is.null(problem$assignmentsLB)) {
    for (i in seq_len(nrow(problem$assignmentsLB))) {
      if (index %in% indices)
        res$assignmentsLB <- rbind(res$assignmentsLB, problem$assignmentsLB[i, ])
      index <- index + 1
    }
  }
  
  if (!is.null(problem$assignmentsUB)) {
    for (i in seq_len(nrow(problem$assignmentsUB))) {
      if (index %in% indices)
        res$assignmentsUB <- rbind(res$assignmentsUB, problem$assignmentsUB[i, ])
      index <- index + 1
    }
  }
  
  if (!is.null(problem$assignmentPairwiseAtLeastComparisons)) {
    for (i in seq_len(nrow(problem$assignmentPairwiseAtLeastComparisons))) {
      if (index %in% indices)
        res$assignmentPairwiseAtLeastComparisons <- rbind(res$assignmentPairwiseAtLeastComparisons,
                                                          problem$assignmentPairwiseAtLeastComparisons[i, ])
      index <- index + 1
    }
  }
  
  if (!is.null(problem$assignmentPairwiseAtMostComparisons)) {
    for (i in seq_len(nrow(problem$assignmentPairwiseAtMostComparisons))) {
      if (index %in% indices)
        res$assignmentPairwiseAtMostComparisons <- rbind(res$assignmentPairwiseAtMostComparisons,
                                                                 problem$assignmentPairwiseAtMostComparisons[i, ])
      index <- index + 1
    }
  }
  
  if (!is.null(problem$minimalClassCardinalities)) {
    for (i in seq_len(nrow(problem$minimalClassCardinalities))) {
      if (index %in% indices)
        res$minimalClassCardinalities <- rbind(res$minimalClassCardinalities,
                                                         problem$minimalClassCardinalities[i, ])
      index <- index + 1
    }
  }
  
  if (!is.null(problem$maximalClassCardinalities)) {
    for (i in seq_len(nrow(problem$maximalClassCardinalities))) {
      if (index %in% indices)
        res$maximalClassCardinalities <- rbind(res$maximalClassCardinalities,
                                                         problem$maximalClassCardinalities[i, ])
      index <- index + 1
    }
  }
  
  return (res)
}

