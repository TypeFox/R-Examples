#' @title Reassign low sample states to close states
#'
#' @description 
#' This function removes small sample states by reassigning points in those 
#' state to nearby states.
#' 
#' This can become necessary when in an iterative algorithm (like 
#' \code{\link{mixed_LICORS}}) the
#' weights start moving away from e.g. state \eqn{j}.  At some point
#' the effective sample size of state \eqn{j} (sum of column \eqn{\mathbf{W}_j}) is so 
#' small that state-conditional estimates (mean, variance, kernel density 
#' estimate, etc.) can not be obtained accurately anymore.  Then it is good to 
#' remove state \eqn{j} and reassign its samples to other (close) states.
#' 
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param min minimum effective sample size to stay in the weight matrix
#' @keywords method manip
#' @export
#' @examples
#' set.seed(10)
#' WW = matrix(c(rexp(1000, 1/10), runif(1000)), ncol =5, byrow=FALSE)
#' WW = normalize(WW)
#' colSums(WW)
#' remove_small_sample_states(WW, 20)

remove_small_sample_states <- function(weight.matrix, min) {
  invisible( normalize(weight.matrix[, colSums(weight.matrix) > min]) )
} 