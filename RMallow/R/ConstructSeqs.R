#' Constructs sequences from Kendall Information matricies.
#' 
#' Sequences in a fully-ordered sequence space have a unique Kendall
#' Information vector associated with them.  This function creates the sequence
#' from the Kendall information vector.
#' 
#' 
#' @param prefs Ordering preference between columns in the data. 1 cooresponds
#' to an increase, 0 to a decrease.
#' @param n.abils Number of columns in the original data set.
#' @return List of fully-ordered sequences, one for each row of prefs.
#' @author Erik Gregory
#' @keywords Sequences
#' @examples
#' ConstructSeqs(matrix(c(1, 1, 1, 0, 0, 0), nrow = 1), 4)
#' # Should output (4, 1, 2, 3)
ConstructSeqs <-
function(prefs, n.abils) {
  R <- list()
  n <- length(prefs)
  singled <- (1 %in% dim(as.matrix(prefs)))
  if (singled) {
    G <- 1:2
    prefs <- rbind(prefs, prefs)
  }
  else {
    G <- 1:ncol(prefs)
  }
  nums <- list()
  for (j in G) {
    nums[[j]] <- 1:n.abils
  }
  tops <- c(0, cumsum(rev(2:n.abils) - 1))
  seqs <- matrix(0, nrow = length(G), ncol = n.abils)
  for (i in 1:(n.abils - 1)) {
    if ((tops[i] + 1) == tops[i + 1]) {
      sumz <- prefs[, ((tops[i] + 1):tops[i + 1])]
    }
    else {
      sumz <- rowSums(prefs[, ((tops[i] + 1):tops[i + 1])])
    }
    for (j in G) {
      abil <- nums[[j]][1 + sumz[j]] 
      if (i == 1) {
        R[[j]] <- abil
      }
      else {
        R[[j]] <- c(R[[j]], abil) 
      }
      nums[[j]] <- nums[[j]][-(1 + sumz[j])]
      if (length(R[[j]]) == (n.abils - 1)) {
        R[[j]] <- c(R[[j]], nums[[j]])
      }
    }
  }
  if(singled) {
    R <- list(R[[1]])
  }
  return(R)
  
}
