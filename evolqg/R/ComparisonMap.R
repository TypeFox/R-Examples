#' Generic Comparison Map functions for creating parallel list methods
#' Internal functions for making eficient comparisons.
#' @param matrix.list list of matrices being compared
#' @param MatrixCompFunc Function used to compare pair of matrices, must output a vector: comparisons and probabilities
#' @param ... Aditional arguments to MatrixCompFunc
#' @param repeat.vector Vector of repeatabilities for correlation correction.
#' @param parallel if TRUE computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return Matrix of comparisons, matrix of probabilities.
#' @import plyr
#' @importFrom reshape2 acast
#' @author Diogo Melo
#' @seealso \code{\link{MantelCor}}, \code{\link{KrzCor}},\code{\link{RandomSkewers}}
ComparisonMap <- function (matrix.list, MatrixCompFunc, ..., repeat.vector = NULL, parallel = FALSE){
  n.matrix <- length(matrix.list)
  if(is.null(names(matrix.list))) {names(matrix.list) <- 1:n.matrix}
  matrix.names <- names (matrix.list)
  CompareToN <- function(n) ldply(matrix.list[(n+1):n.matrix],
                                  function(x) {MatrixCompFunc(x, matrix.list[[n]], ...)[1:2]},
                                  .parallel = parallel)
  comparisons <- adply(1:(n.matrix-1), 1,  CompareToN, .parallel = parallel)
  corrs <- suppressMessages(acast(comparisons[-4], X1~.id)[,matrix.names[-1]])
  probs <- suppressMessages(acast(comparisons[-3], X1~.id)[,matrix.names[-1]])
  probabilities <- array (0, c(n.matrix, n.matrix))
  correlations <- probabilities
  probabilities[upper.tri(probabilities)] <- probs[upper.tri(probs, diag=T)]
  correlations[upper.tri(correlations)] <- corrs[upper.tri(probs, diag=T)]
  if (!is.null (repeat.vector)) {
    repeat.matrix <- sqrt(outer(repeat.vector, repeat.vector))
    correlations[lower.tri(correlations)] <- t(correlations/repeat.matrix)[lower.tri(correlations)]
    diag (correlations) <- repeat.vector
  }
  rownames (correlations) <- matrix.names
  colnames (correlations) <- matrix.names
  dimnames (probabilities) <- dimnames (correlations)
  output <- list ('correlations' = t(correlations), 'probabilities' = t(probabilities))
  return (output)
}
