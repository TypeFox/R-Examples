#' Find Indices for Non-identifiable Elements of State Matrix.
#'
#' \code{find_wildcards} finds indices for non-identifiable elements of state
#' matrix.
#'
#' This is a helper function for \code{\link{degeneracy_check}}.
#'
#' @param state_mat Numeric matrix with rows as states and columns as
#'   predictors.
#' @param action_vec Numeric vector indicating what action to take for each
#'   state.
#' @param cols Numeric vector same length as number of columns of the
#'   state matrix (\code{state_mat}) with the action that each column of the
#'   state matrix corresponds to the decision model taking in the previous
#'   period. This is only relevant when the predictor variables of the FSM are
#'   lagged outcomes that include the previous actions taken by that decision
#'   model.
#'
#' @return Returns a list of indices (tuples specifying row and column of a
#'   matrix).
#'
#' @examples
#' tft_state <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2), 2, 4)
#' tft_action <- matrix(c(1, 2))
#' find_wildcards(tft_state, tft_action,  c(1, 2, 1, 2))
#'
#' @export

find_wildcards <- function(state_mat, action_vec, cols){
  # cols is a vector that says what action each col of state_mat corresponds to
  if (ncol(state_mat)!=length(cols)){
    stop("Number of cols of state matrix
                     does not equal length of cols argument.")
  }
  if (nrow(state_mat)!=length(action_vec)){
    stop("Number of rows of state matrix does not equal
                     the length of the action vector.")
  }

  counter <- 1
  indices <- as.list(rep(NA, length(as.vector(state_mat))))
  for (i in seq(length(action_vec))){
    for (j in seq(length(cols))){
      if (action_vec[i] != cols[j]){
        indices[[counter]] <- c(i, j)
        counter <- counter + 1
      } else indices[[counter]] <- NULL
    }
  }
  indices
}


# cols <- c(1, 2, 1, 2)

#' Determines if State Matrix is Degenerate for Given Data Set.
#'
#' \code{degeneracy_check} finds indices for non-identifiable elements of state
#' matrix and then flips values for those elements and checks changes in
#' resulting fitness.
#'
#' \code{degeneracy_check} finds indices for non-identifiable elements of state
#' matrix and then flips values for those elements and checks changes in
#' resulting fitness. Being in state/row k (e.g. 2) corresponds to taking action
#' j (e.g. D). For row k, all entries in the matrix that corresponds to taking
#' action j last period (e.g. columns 2 and 4 for D) are identifiable; however,
#' columns that correspond to not taking action j last period (e.g. columns 1
#' and 3 for D) for the row $k$ that corresponds to taking action j are not
#' identifiable for a deterministic play of the strategy. For all elements of
#' the matrix that are not identifiable, the value of the element can be any
#' integer in the inclusive range of the number of rows of the matrix (e.g. 1 or
#' 2). With empirical data, where the probability that a single deterministic
#' model generated the data is effectively zero, it is useful to find every
#' entry in the matrix that would be unidentifiable if the strategy were played
#' deterministically and then for each element flip it to its opposite value and
#' test for any change in fitness of the strategy on the data. This function
#' implements this idea. If there is no change, a sparse matrix is returned
#' where the the elements in that matrix with a 0 are unidentifiable because
#' their value makes no difference to the fit of the strategy to the provided
#' data. If, for each element in the matrix, switching its value led to a
#' decrease in fitness the following message is displayed, ``Your strategy is a
#' deterministic approximation of a stochastic process and all of the elements
#' of the state matrix can be identified.'' If the model is fine, then
#' \code{sparse_state_mat} and \code{corrected_state_mat} should be equal to
#' \code{state_mat}.
#'
#' @param state_mat Numeric matrix with rows as states and columns as
#'   predictors.
#' @param action_vec Numeric vector indicating what action to take for each
#'   state.
#'  @param cols Optional numeric vector same length as number of columns of the
#'   state matrix (\code{state_mat}) with the action that each column of the
#'   state matrix corresponds to the decision model taking in the previous
#'   period. This is only relevant when the predictor variables of the FSM are
#'   lagged outcomes that include the previous actions taken by that decision
#'   model.
#' @param data Numeric matrix that has first col period and rest of cols are
#'   predictors.
#' @param outcome Numeric vector same length as the number of rows as data.
#'
#' @return Returns a list of with sparse and corrected state matrix.  If the
#'   model is fine, then \code{sparse_state_mat} and \code{corrected_state_mat}
#'   should be equal to \code{state_mat}.
#'
#' @export

degeneracy_check <- function(state_mat, action_vec, cols, data, outcome){
  indices <- find_wildcards(state_mat, action_vec, cols)

  results1 <- fitnessCPP(action_vec, state_mat, data)
  if (anyNA(results1) | length(results1)==0){
    stop("Results from first fitness evaluation have missing values.")
  }
  results1 <- sum(ifelse( results1 == outcome , 1 , 0)) / length(results1)

  sparse_state_mat <- state_mat
  corrected_state_mat <- state_mat
  # sparse will be corrected and have zeros added for non-identifiable
  # corrected will only be corrected, this way it can be directly used for varImp()

  dif <- rep(NA, length(indices)) # how many are not identifiable

  for(i in seq(length(dif))){

    state_mat_flipped <- state_mat
    state_mat_flipped[indices[[i]][1],
                      indices[[i]][2]] <- ifelse(state_mat[indices[[i]][1],
                                                           indices[[i]][2]]==1, 2, 1)

    results2 <- fitnessCPP(action_vec, state_mat_flipped, data)
    if (anyNA(results2) | length(results2)==0){
      stop("Results from subsequent fitness evaluation have missing values.")
    }
    results2 <- sum(ifelse( results2 == outcome , 1 , 0)) / length(results2)

    dif[i] <- (results2 - results1) / results1

    # for any equal zero, use 0 as the identifier of non-identifiable
    sparse_state_mat[indices[[i]][1],
                     indices[[i]][2]] <- ifelse(dif[i]==0, 0, state_mat[indices[[i]][1],
                                                                        indices[[i]][2]])
    # for any greater than zero, give it opposite value that improved it
    sparse_state_mat[indices[[i]][1],
                     indices[[i]][2]] <- ifelse(dif[i] > 0, state_mat_flipped[indices[[i]][1],
                                                                              indices[[i]][2]], state_mat[indices[[i]][1],
                                                                                                          indices[[i]][2]])

    # for any greater than zero, give it opposite value that improved it
    corrected_state_mat[indices[[i]][1],
                        indices[[i]][2]] <- ifelse(dif[i] > 0, state_mat_flipped[indices[[i]][1],
                                                                                 indices[[i]][2]], state_mat[indices[[i]][1],
                                                                                                             indices[[i]][2]])
  }
  # if the model is fine, then sparse_state_mat and corrected_state_mat should be equal to state_mat
  list(dif = dif, sparse_state_mat = sparse_state_mat, corrected_state_mat = corrected_state_mat)
}

# if(any(dif > 0)) cat("You have not found an optimal strategy.\n")
# if(any(dif < 0)) cat("You have randomness in your strategy, which may or may not be optimal.\n")
