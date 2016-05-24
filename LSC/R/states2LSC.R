#' @title Estimate local statistical complexity (LSC) from states
#'
#' @description 
#' Converts states (either as a size \eqn{N} array of labels or an 
#' \eqn{N \times K} weight matrix) into \eqn{N} local statistical complexities 
#' (LSC) per state.
#' 
#' If \code{states} is a matrix - representing the state space in the 
#' dimensions of the original data - then the LSC output will be formatted 
#' automatically to an array of the same shape and dimension.
#' 
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param states array of size \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}.
#' @param base logarithm base for complexity (entropy). Default: \code{base = 2}
#'  (thus 'bits').
#' @param type estimation type for the probabilities: \code{"MLE"}
#' @keywords manip
#' @export
#' @seealso \code{\link{states2probs}}
#'
states2LSC <- function(weight.matrix = NULL, states = NULL,  
                       base = 2, type = c("MLE")) {
  
  lsc <- -log(states2probs(states = states, weight.matrix = weight.matrix, type = type), 
              base)
  invisible(lsc) 
}
