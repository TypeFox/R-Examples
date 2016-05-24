#' Decodes Action Vector
#'
#' \code{decode_action_vec} decodes action vector.
#'
#' This function takes a solution string of binary values in Gray
#' representation, transforms it to a decimal representation, then puts it in
#' matrix form with the correct sized matrices, given the specified numbers of
#' states, inputs, and actions.
#'
#' @param string Numeric (integer) vector of only 1's and 0's.
#' @param states Numeric vector with the number of states, which is the number
#'   of rows.
#' @param inputs Numeric vector length one, with the number of columns.
#' @param actions Numeric vector with the number of actions. Actions (and
#'   states) determine how many binary elements we need to represent an element
#'   of the action (or state) matrix.
#'
#' @return Returns numeric (integer) vector.
#'
#' @export

decode_action_vec <- function(string, states, inputs, actions){
  # "states" is the number of rows;  "inputs" is the number of columns
  # "actions" (and "states") specifies how many binary elements we need to represent
  # an element of the action (or state) matrix.

  poss.state.values <- (1:states) - 1
  b1 <- GA::decimal2binary(max(poss.state.values))
  l1 <- length(b1) #how many binary elements to represent one element of state matrix

  poss.action.values <- (1:actions) - 1
  b2 <- GA::decimal2binary(max(poss.action.values))
  l2 <- length(b2) #how many binary elements to represent one element of action matrix

  #        if (length(string) != states * inputs * l1 + states * l2) {
  #                cat('string = ', as.character(string), ', length = ', length(string), '\n')
  #                cat('states = ', states, ', inputs = ', inputs, ', actions = ', actions, '\n')
  #                cat('l1 = ', l1, ', l2 = ', l2, ', target length = ', states * inputs * l1 + states * l2, '\n')
  #        }
  stopifnot(length(string) == (states*inputs*l1 + states*l2))
  # this is states*l2 because action.vec has length == nrows(state.matrix)

  string2vec <- function(string, len){
    vec_len <- as.integer(length(string) / len)
    string <- GA::gray2binary(string)

    vec <- integer(vec_len)

    for (i in 1:vec_len){
      idx <- (i-1) * len
      vec[i] <- GA::binary2decimal(string[idx + 1:len])
    }
    as.vector(vec) + 1
  }

  action.string <- string[1:(states*l2)]
  action.vec <- string2vec(action.string, l2)
  # When you represent decimals in binary form with 2 binaries
  # this can represent 0, 1, 2, or 3, so we use to change all 3s to 2s because here
  # we needed only 1 and 2 for the elements of both our action and state matrices.
  #        action.vec[action.vec==3] <- 2
  #        action.vec[action.vec==0] <- 1
  as.integer(action.vec)
}
