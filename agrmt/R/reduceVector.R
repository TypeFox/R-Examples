reduceVector <-
function(X) {
  # Helper function to remove all zeros, and then remove repeated values
  # Argument: vector
  X <- X[!X==0]             # remove all 0 from X
  n <- length(X)
  Z <- X[1]                 # first item of X is set in Z
  for (i in 2:n) {          # for all subsequent items in Z
    if (!X[i] == X[i-1]) Z <- c(Z,X[i]) # if items is not the same, add it to Z
    }
  return(Z)
  }
