#' Builds Bitstring
#'
#' \code{build_bitstring} creates a bitstring from an action vector, state
#' matrix, and number of actions.
#'
#'
#' @param action_vec Numeric vector indicating what action to take for each
#'   state.
#' @param state_mat Numeric matrix with rows as states and columns as
#'   predictors.
#' @param actions Numeric vector length one with the number of actions.
#'
#' @return Returns numeric vector bitstring.
#' @export

build_bitstring <- function(action_vec, state_mat, actions) {
  states <- nrow(state_mat)
  inputs <- ncol(state_mat)

  if (length(action_vec) != states)
    stop("Error: action_vec does not match number of states")

  poss.state.values <- (1:states) - 1
  b1 <- GA::decimal2binary(max(poss.state.values))
  l1 <- length(b1) #how many binary elements to represent one element of state matrix

  poss.action.values <- (1:actions) - 1
  b2 <- GA::decimal2binary(max(poss.action.values))
  l2 <- length(b2) #how many binary elements to represent one element of action matrix

  bitstring <- integer(states*inputs*l1 + states*l2)
  bitstring[1:(l2 * states)] <- GA::binary2gray(unlist(lapply(action_vec,
                                                              function(x) GA::decimal2binary(x-1,l2))))
  bitstring[(1 + l2*states):length(bitstring)] <-
          GA::binary2gray(unlist(lapply(as.vector(state_mat),
                                        function(x) GA::decimal2binary(x-1,l1))))

  bitstring
}
