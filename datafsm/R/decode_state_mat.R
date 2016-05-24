#' Decodes State Matrix
#'
#' \code{decode_state_mat} decodes state matrix.
#'
#' This function takes a solution string of binary values in Gray
#' representation, transforms it to a decimal representation, then puts it in
#' matrix form with the correct sized matrices, given the specified numbers of
#' states, inputs, and actions.
#'
#' @param string Numeric vector.
#' @param states Numeric vector with the number of states, which is the number
#'   of rows.
#' @param inputs Numeric vector length one, with the number of columns.
#' @param actions Numeric vector with the number of actions. Actions (and
#'   states) determine how many binary elements we need to represent an element
#'   of the action (or state) matrix.
#'
#' @return Returns numeric (integer) matrix.
#'
#' @export

decode_state_mat <- function(string, states, inputs, actions){
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

  string2matrix <- function(string, len){
    vec_len <- as.integer(length(string) / len)
    string <- GA::gray2binary(string)

    vec <- integer(vec_len)


    for (i in 1:vec_len){
      idx <- (i-1) * len
      vec[i] <- GA::binary2decimal(string[idx + 1:len])
    }
    vec <- as.vector(vec) + 1

    matrix(vec,
           nrow = states, ncol = inputs, byrow = FALSE)
  }

  state.string <- string[(states*l2+1):length(string)]

  state.matrix <- string2matrix(state.string, l1)
  state.matrix
}
