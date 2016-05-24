#' @title Converts label vector to 0/1 mixture weight matrix 
#' 
#' @description 
#' Converts unique cluster assignment stored in a length \eqn{N} label vector 
#' into a \eqn{N \times K} Boolean matrix of mixture weights.
#' 
#' @param states a vector of length \eqn{N} with the state labels
#' @param num.states.total total number of states. If \code{NULL}, then the 
#' maximum of \code{states} is chosen
#' @keywords manip array
#' @export
#' @seealso \code{\link{weight_matrix2states}}
#' @examples
#' ss = sample.int(5, 10, replace = TRUE)
#' WW = states2weight_matrix(ss)
#' 
#' image2(WW, col = "RdBu", xlab = "States", ylab = "Samples", axes = FALSE)
#'       

states2weight_matrix <- function(states, num.states.total = NULL){
  
  if (is.null(num.states.total)){
    num.states.total <- max(states)
  }
  weight.matrix <- Matrix(0, nrow = length(states), ncol = num.states.total,
                          sparse = TRUE)
  weight.matrix[cbind(seq_along(states), states)] <- 1
  invisible(weight.matrix)
}

