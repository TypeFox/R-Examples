#' Transition Matrix Estimation
#'
#' This function allows you to estimate transition matrices (probabilities or counts) for up-to-fifth-order discrete Markov chains. For n-order Markov chains with n greater than 1, you can access the estimated transition matrices through nested lists.
#' @param sequence A vector of integers representing the sequence. The sequence must be in the form of 1,2,3,...
#' @param order Integer from 1 to 5. Order of the Markov chain. By default is set to 1.
#' @param probs Logical. If TRUE probability matrices are returned, otherwise count matrices are returned. By default is set to TRUE.
#' @keywords markov chain transition matrix
#' @export TransMatrix
#' @examples
#' seq <- sample(c(1,2,3,4), size = 1000, replace = TRUE)
#' TransMatrix(seq, order = 1, probs = TRUE)
#' TransMatrix(seq, order = 2, probs = FALSE)
#' mc <- TransMatrix(seq, order = 4, probs = TRUE)
#' mc[[1]][[2]][[3]] # through nested lists you can access to the estimated transition matrices

TransMatrix <- function(sequence, order = 1, probs = TRUE) {

  # ORDER 1
  if (order == 1) {
    states <- sort(unique(sequence))
    n.states <- length(states)
    trans.mat <- matrix(0, nrow = n.states, ncol = n.states)
    colnames(trans.mat) <- rownames(trans.mat) <- states

    for (i in 1:n.states) {
      for (j in 1:n.states) {
        idx <- which(sequence == states[i])
        idx <- idx[idx < length(sequence)]
        trans.mat[i,j] <- length(which(sequence[idx+1] == states[j]))
      }
    }

    # RETURN
    if (probs)
      return(trans.mat/rowSums(trans.mat))
    else
      return(trans.mat)
  }

  # ORDER 2
  if (order == 2) {
    states <- sort(unique(sequence))
    n.states <- length(states)
    trans.mat <- rep(list( matrix(0, nrow = n.states, ncol = n.states) ), n.states)

    for (k in 1:n.states) {
      colnames(trans.mat[[k]]) <- rownames(trans.mat[[k]]) <- states
      for (i in 1:n.states) {
        for (j in 1:n.states) {
          idx <- which(sequence == states[k])
          idx <- idx[idx < (length(sequence)-1)]
          trans.mat[[k]][i,j] <- length(which(sequence[idx+1] == states[i] &
                                              sequence[idx+2] == states[j]))
        }
      }
    }

    # RETURN
    if (probs) {
      for (k in 1:n.states) {
        trans.mat[[k]] <- trans.mat[[k]] / rowSums(trans.mat[[k]])
      }
      return(trans.mat)
    } else {
        return(trans.mat)
    }
  }

  # ORDER 3
  if (order == 3) {
    states <- sort(unique(sequence))
    n.states <- length(states)
    trans.mat <- list()

    for (h in 1:n.states) {
      trans.mat[[h]] <- rep(list( matrix(0, nrow = n.states, ncol = n.states) ), n.states)
      for (k in 1:n.states) {
        colnames(trans.mat[[h]][[k]]) <- rownames(trans.mat[[h]][[k]]) <- states
        for (i in 1:n.states) {
          for (j in 1:n.states) {
            idx <- which(sequence == states[k])
            idx <- idx[idx < (length(sequence)-2)]
            trans.mat[[h]][[k]][i,j] <- length(which(sequence[idx+1] == states[k] &
                                                sequence[idx+2] == states[i] &
                                                sequence[idx+3] == states[j]))
          }
        }
      }
    }

    # RETURN
    if (probs) {
      for (h in 1:n.states) {
        for (k in 1:n.states) {
          trans.mat[[h]][[k]] <- trans.mat[[h]][[k]] / rowSums(trans.mat[[h]][[k]])
        }
      }
      return(trans.mat)
    } else {
      return(trans.mat)
    }
  }

  # ORDER 4
  if (order == 4) {
    states <- sort(unique(sequence))
    n.states <- length(states)
    trans.mat <- list()

    for (v in 1:n.states) {
      trans.mat[[v]] <- list()
      for (h in 1:n.states) {
        trans.mat[[v]][[h]] <- rep(list( matrix(0, nrow = n.states, ncol = n.states) ),
                                   n.states)
          for (k in 1:n.states) {
            colnames(trans.mat[[v]][[h]][[k]]) <- rownames(trans.mat[[v]][[h]][[k]]) <- states
            for (i in 1:n.states) {
              for (j in 1:n.states) {
                idx <- which(sequence == states[k])
                idx <- idx[idx < (length(sequence)-3)]
                trans.mat[[v]][[h]][[k]][i,j] <- length(which(sequence[idx+1] == states[h] &
                                                       sequence[idx+2] == states[k] &
                                                       sequence[idx+3] == states[i] &
                                                       sequence[idx+4] == states[j]))
              }
            }
          }
        }
      }

    # RETURN
    if (probs) {
      for (v in 1:n.states) {
        for (h in 1:n.states) {
          for (k in 1:n.states) {
            trans.mat[[v]][[h]][[k]] <- trans.mat[[v]][[h]][[k]] / rowSums(trans.mat[[v]][[h]][[k]])
          }
        }
      }
      return(trans.mat)
    } else {
      return(trans.mat)
    }
  }

  # ORDER 5
  if (order == 5) {
    states <- sort(unique(sequence))
    n.states <- length(states)
    trans.mat <- list()

    for (z in 1:n.states) {
      trans.mat[[z]] <- list()
      for (v in 1:n.states) {
        trans.mat[[z]][[v]] <- list()
        for (h in 1:n.states) {
          trans.mat[[z]][[v]][[h]] <- rep(list( matrix(0, nrow = n.states, ncol = n.states) ),
                                   n.states)
          for (k in 1:n.states) {
            colnames(trans.mat[[z]][[v]][[h]][[k]]) <- rownames(trans.mat[[z]][[v]][[h]][[k]]) <- states
            for (i in 1:n.states) {
              for (j in 1:n.states) {
                idx <- which(sequence == states[k])
                idx <- idx[idx < (length(sequence)-5)]
                trans.mat[[z]][[v]][[h]][[k]][i,j] <- length(which(sequence[idx+1] == states[v] &
                                                              sequence[idx+2] == states[h] &
                                                              sequence[idx+3] == states[k] &
                                                              sequence[idx+4] == states[i] &
                                                              sequence[idx+5] == states[j]))
              }
            }
          }
        }
      }
    }

    # RETURN
    if (probs) {
      for (z in 1:n.states) {
        for (v in 1:n.states) {
          for (h in 1:n.states) {
            for (k in 1:n.states) {
              trans.mat[[z]][[v]][[h]][[k]] <- trans.mat[[z]][[v]][[h]][[k]] / rowSums(trans.mat[[z]][[v]][[h]][[k]])
            }
          }
        }
      }
      return(trans.mat)
    } else {
      return(trans.mat)
    }
  }
}

