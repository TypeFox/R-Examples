collapse <-
  function(D, pos=FALSE) {
    # Calculates a frequency vector F from a population vector D
    # Arguments:  D  = population vector
    #            pos = positions of the categories (necessary if there are categories with 0 observations)
    #                 if nothing i provided, we assume observations in all categories
    #                 e.g. if position 2 is not present in the population, no 0 will be added
    if(is.vector(D) == FALSE) {
      warning("Warning: Expected a vector, or vector is empty.")
      return(NA)}
    T <- as.data.frame(table(D))  # table(D) to count frequencies
    F <- T[,2]                    # [,2] chooses the values
    l <- length(pos)              # number of categories specified
    if (!(pos[1] == FALSE & pos[l] == FALSE)) {   # number of positions is provided:
      # if first value of pos is not zero, and last value is not zero
      F2 <- F           # source of values
      k <- length(pos)  # number of positions
      F <- rep(0,k)     # empty with length(k)
      V <- unique(D)    # chooses the categories (positions) of the values in F2
      V <- V[!is.na(V)] # remove NA from unique values
      V <- sort(V)      # keeping order that is messed up by the preceding line
      m <- length(V)    # number of non-zero values to be filled into F
      for(i in 1:m) {
        for(j in 1:k) {
          if (pos[j] == V[i]) F[j] <- F2[i] # fill in frequencies at the correct positions
        }
      }
    }
    return(F)
  }
