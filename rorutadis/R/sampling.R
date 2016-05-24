#' Stochastic results
#'
#' The function calculates stochastic results for alternative assignments,
#' assignment-based preference relation and class cardinalities.
#' The results are computed by sampling the space of compatible models.
#' 
#' @param problem A problem to consider.
#' @param nrSamples Number of samples. Use more for better quality of results.
#' @return List with the following named elements:
#' \itemize{
#' \item \emph{assignments} - \emph{n} x \emph{p} matrix, where \emph{n} is
#' the number of alternatives and \emph{p} is number of classes; each element
#' \code{[i, j]} contains the rate of samples, for which alternative \emph{a_i}
#' was assigned to class \emph{C_j}.
#' The exact result can be calculated with function \link{calculateAssignments}.
#' \item \emph{preferenceRelation} - \emph{n} x \emph{n} matrix, where \emph{n} is
#' the number of alternatives; each element \code{[i, j]} contains the rate
#' of samples, for which alternative \emph{a_i} was assigned to class at least
#' as good as class of \emph{a_j}.
#' The exact result can be calculated with function \link{compareAssignments}.
#' \item \emph{classCardinalities} - \emph{p} x \emph{(n + 1)} matrix, where \emph{n}
#' is the number of alternatives and \emph{p} is number of classes; each element
#' \code{[i, j]} contains the rate of samples, for which \emph{j-1} alternatives
#' were assigned to class \emph{C_i}. \strong{Note!} first column corresponds to
#' \strong{0} elements.
#' The exact result can be calculated with function \link{calculateExtremeClassCardinalities}.
#' }
#' @seealso
#' \code{\link{buildProblem}}
#' \code{\link{calculateAssignments}}
#' \code{\link{compareAssignments}}
#' \code{\link{calculateExtremeClassCardinalities}}
#' @examples
#' perf <- matrix(c(2,1,1,2), 2)
#' problem <- buildProblem(perf, 2, FALSE, c('g', 'g'), c(0, 0))
#' 
#' calculateStochasticResults(problem, 1000)
#' @import hitandrun
#' @export
calculateStochasticResults <- function(problem, nrSamples = 100) {
  stopifnot(nrSamples > 0)
  
  model <- buildModel(problem, TRUE)
  
  if (!isModelConsistent(model)) {
    stop("Model infeasible.")
  }
  
  isSpaceContinuous <- all(model$constraints$types == "C")
  
  if (!isSpaceContinuous) {
    warning("Calculating stochastic results for problems with assignement-based pairwise 
            comparisons and/or desired class cardinalities may take a lot of time 
            due to discontinuous character of sampling space.", immediate. = TRUE)
  }
  
  nrAlternatives <- nrow(problem$perf)
  nrCriteria <- ncol(problem$perf)
  nrClasses <- problem$nrClasses
  
  result <- list(assignments = NULL, preferenceRelation = NULL, classCardinalities = NULL)
  
  result$assignments <- matrix(data = 0, nrow = nrAlternatives, ncol = nrClasses)
  result$preferenceRelation <- matrix(data = 0, nrow = nrAlternatives, ncol = nrAlternatives)
  result$classCardinalities <- matrix(data = 0, nrow = nrClasses, ncol = nrAlternatives + 1)
  
  if (isSpaceContinuous) {
    model <- buildModel(problem, FALSE)
  } else {
    problemWithOnlyAssignmentsAsPreferenceInformation <- problem
    problemWithOnlyAssignmentsAsPreferenceInformation$assignmentPairwiseAtLeastComparisons <- NULL
    problemWithOnlyAssignmentsAsPreferenceInformation$assignmentPairwiseAtMostComparisons <- NULL
    problemWithOnlyAssignmentsAsPreferenceInformation$minimalClassCardinalities <- NULL
    problemWithOnlyAssignmentsAsPreferenceInformation$maximalClassCardinalities <- NULL
    
    problemWithOnlyAssignmentsAsPreferenceInformation
    model <- buildModel(problemWithOnlyAssignmentsAsPreferenceInformation, FALSE)
  }
  
  constraints <- model$constraints
  constraints$dir[which(constraints$dir == "==")] <- "="
  geq <- which(constraints$dir == ">=")
  
  for (i in geq) {
    constraints$rhs[i] <- -1 * constraints$rhs[i]
    constraints$lhs[i, ] <- -1 * constraints$lhs[i, ]
  }
  
  constraints$dir[geq] <- "<="
  names(constraints)[1] <- "constr"
  constraints[[4]] <- NULL
  
  state <- har.init(constraints, thin.fn = function(n) { ceiling(log(n + 1)/4 * n^3) },
                    thin = NULL, x0.randomize = FALSE, x0.method = "slacklp", x0 = NULL)
  
  producedSamples <- 0
  allGeneratedSamples <- 0
  
  while (producedSamples < nrSamples) {
    harSample <- har.run(state, n.samples = 1)
    state <- harSample$state
    sample <- harSample$samples[1, ]
    
    thresholds <- getThresholdsFromF(model, sample)
    assignments <- getAssignmentsFromF(model, sample, thresholds)
    
    allGeneratedSamples <- allGeneratedSamples + 1
    
    if (isSpaceContinuous || isAssignmentsValid(problem, assignments)) { # isAssignmentsValid makes rejection sampling for discontinuous spaces
      for (i in seq_len(length(assignments))) {
        result$assignments[i, assignments[i]] <- result$assignments[i, assignments[i]] + 1
      }
      
      for (i in seq_len(nrAlternatives)) {
        for (j in i:nrAlternatives) {
          if (i == j) {
            result$preferenceRelation[i, j] <- result$preferenceRelation[i, j] + 1
          }
          else if (assignments[i] == assignments[j]) {
            result$preferenceRelation[i, j] <- result$preferenceRelation[i, j] + 1
            result$preferenceRelation[j, i] <- result$preferenceRelation[j, i] + 1
          }
          else if (assignments[i] > assignments[j]) {
            result$preferenceRelation[i, j] <- result$preferenceRelation[i, j] + 1
          }
          else {
            result$preferenceRelation[j, i] <- result$preferenceRelation[j, i] + 1
          }
        }
      }
      
      cardinalities <- sapply(seq_len(nrClasses), function(h) { length(which(assignments == h)) })
      
      for (h in seq_len(nrClasses)) {
        result$classCardinalities[h, cardinalities[h] + 1] <- result$classCardinalities[h, cardinalities[h] + 1] + 1
      }
      
      producedSamples <- producedSamples + 1
    }
  }
  
  result$assignments <- result$assignments / producedSamples
  result$preferenceRelation <- result$preferenceRelation / producedSamples
  result$classCardinalities <- result$classCardinalities / producedSamples
  result$acceptanceRateOfSecondPhaseSampling <- producedSamples / allGeneratedSamples
    
  return (result)
}

isAssignmentsValid <- function(problem, assignmentsToCheck) {
  for (i in seq(nrow(problem$assignmentsLB))) {
    a <- problem$assignmentsLB[i, 1]
    
    if (assignmentsToCheck[a] < problem$assignmentsLB[i, 2])
      return (FALSE)
  }
  
  for (i in seq(nrow(problem$assignmentsUB))) {
    a <- problem$assignmentsUB[i, 1]
    
    if (assignmentsToCheck[a] > problem$assignmentsUB[i, 2])
      return (FALSE)
  }
  
  for (i in seq(nrow(problem$assignmentPairwiseAtLeastComparisons))) {
    a <- problem$assignmentPairwiseAtLeastComparisons[i, 1]
    b <- problem$assignmentPairwiseAtLeastComparisons[i, 2]
    
    if (assignmentsToCheck[a] - assignmentsToCheck[b] < 
          problem$assignmentPairwiseAtLeastComparisons[i, 3])
      return (FALSE)
  }
  
  for (i in seq(nrow(problem$assignmentPairwiseAtMostComparisons))) {
    a <- problem$assignmentPairwiseAtMostComparisons[i, 1]
    b <- problem$assignmentPairwiseAtMostComparisons[i, 2]
    
    if (assignmentsToCheck[a] - assignmentsToCheck[b] > 
          problem$assignmentPairwiseAtMostComparisons[i, 3])
      return (FALSE)
  }
  
  for (i in seq(nrow(problem$minimalClassCardinalities))) {
    h <- problem$minimalClassCardinalities[i, 1]
    
    if (length(which(assignmentsToCheck == h)) < problem$minimalClassCardinalities[i, 2])
      return (FALSE)
  }
  
  for (i in seq(nrow(problem$maximalClassCardinalities))) {
    h <- problem$maximalClassCardinalities[i, 1]
    
    if (length(which(assignmentsToCheck == h)) > problem$maximalClassCardinalities[i, 2])
      return (FALSE)
  }
  
  return (TRUE)
}
