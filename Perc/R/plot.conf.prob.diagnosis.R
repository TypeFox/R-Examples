# plotProbDiagnosis: binary one, allow users to specify the value, ordered by sim.rank
#' Diagnosis Plot
#' \code{plotProbDiagnosis} generate heat map for dominance probability matrix
#' @param prob.mat dominance probability matrix
#' @param cutoff a numeric value between 0.5 to 1. 
#' A value that is equal or greater than the cutoff is considered of high certainty.
#' @param ... Further argument may be supplied and processed by \code{levelplot}.
#' 
#' @seealso \code{\link{plotConfmat}}

plotProbDiagnosis <- function(prob.mat, cutoff = 0.75, ...) {
  if (!(all(prob.mat >= 0) & all(prob.mat <= 1)))
    stop("Use the Dominance Probability matrix; if you want to plot any matrix, use plot.conf.mat")
  if (!(cutoff >= 0.5 & cutoff <= 1))
    stop("A cutoff value is no smaller than 0.5 and no greater than 1.")
  if (min(prob.mat) < 0.5) probMatrix <- valueConverter(prob.mat)
  matrix <- myrecode(
    probMatrix, 
    list(
      probMatrix[probMatrix < cutoff], 
      probMatrix[probMatrix >=  cutoff]
      ),
    c(0, 1)
    )
  
  plotConfmat(matrix, ordering = NA, labels = FALSE, col.regions = c("white", "black"))
}