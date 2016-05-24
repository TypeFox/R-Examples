#' @title Convert states to vector of probabilities of any given state
#'
#' @description 
#' Converts states (either as a size \eqn{N} array of labels or an 
#' \eqn{N \times K} weight matrix) into \eqn{N} probabilities per state.
#' 
#' If \code{states} is a matrix - representing the state space in the 
#' dimensions of the original data - then the probabilities will be formatted 
#' automatically to an array of the same shape and dimension.
#' 
#' @param states array of size \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param type estimation type for the probabilities: \code{c("MLE")}
#' @keywords manip
#' @export
#' @seealso \code{\link[LICORS]{weight_matrix2states}}
#' @examples
#' 
#' state.sim <- sample.int(5, 100, replace = TRUE)
#' prob.state <- states2probs(states = state.sim)
#' layout(matrix(1:2, ncol = 2))
#' plot(state.sim, xlab = "", ylab = "state")
#' plot(prob.state, xlab = "", ylab = "probability")
#' plot(state.sim, prob.state, xlab = "state", ylab = "probability", type = 'h')
#' 

states2probs <- function(weight.matrix = NULL, states = NULL, type = "MLE") {
  
  if (is.null(states) & is.null(weight.matrix)) {
    stop("You must provide either state vector or the weight matrix.")
  }
  type <- match.arg(type)
  
  if (!is.null(weight.matrix)) {
    marginal.state.probs <- estimate_state_probs(weight.matrix)
    point.probs <- rowSums(weight.matrix %*% Diagonal(x = marginal.state.probs))
      # rowSums(sweep(weight.matrix, 2, marginal.state.probs, "*"))
  } else {
    point.probs <- estimate_state_probs(states = states)[states]
    names(point.probs) <- NULL
  }
  
  if (!is.null(states) & !is.null(dim(states))) {
    invisible(array(point.probs, dim = dim(states)))
  } else {
    invisible(point.probs)
  }
}