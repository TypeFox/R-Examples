agreementError <- function(V, n=500, e=0.01, pos=FALSE) {
  # arguments: V = vector of values (*not* collapsed)
  #            n = number of samples drawn
  #            e = proportion of samples that are errors
  #          pos = vector of possible positions
  #                if FALSE, values in V are set as possible values

  # input validation: is it a vector?
  if(is.vector(V) == FALSE) {
    warning("Warning: Expected a vector, or vector is empty.")
    return(NA)
  }
  
  # input validation: is the position argument provided?
  # if position argument is not provided, then use the values that occur in the vector as positions
  l <- length(pos)
  if (pos[1] == FALSE & pos[l] == FALSE) pos <- as.numeric(as.data.frame(table(V))[,1])

  # input validation: is the pos argument realistic?
  if (length(pos) < 3) {
    warning("Too few positions to calculate agreement: pos argument ignored.")
    return(NA)
  }

  # define samples to draw from, based on input
  nActual <- round(n*(1-e),0) # how many samples to draw *without* error
  nError  <- n - nActual      # how many samples to draw *with* error
  lenPos  <- length(pos)      # how many positions? (provided or assumed)

  # [1] sample nActual values from V (with replacement)
  sActual <- sample(V, nActual, replace=TRUE)

  # [2] sample nError value from V (with replacement)
  # but add or remove 1 (whilst checking if result is within bounds)
  sError <- NULL              # empty to start with
  for(i in 1:nError){
    s <- sample(V, 1, replace=TRUE) # pick a random position
    # add error:
    m <- match(s, pos)        # at which position of pos is the sampled value
    e <- sample(c(-1,0,1),1)  # simulated error
    x <- pos[m+e]             # matched position +/- 1 or exact
      if(m==1) {              # if furthest to left
      e1 <- sample(c(0,1),1)  # censored: cannot become more negative
      x <- pos[m+e1]          # matched position +1 or exact
    }
    if(m==lenPos) {           # if furthest to right
      e2 <- sample(c(-1,0),1) # censored: cannot become more positive
      x <- pos[m+e2]          # matched position -1 or exact
    }
    sError <- c(sError,x)     # add modified (with error) to sample of errors
  }

  # [3] combine samples, and calculate agreement
  sCombined <- c(sActual, sError)
  a <- agreement(collapse(sCombined, pos=pos))
  # now I simply need to run error.agreement() a few (hundred) times
  # and can run the summary statistics on the result (mean, median, SD etc.)
  # z <- replicate(1000, error.agreement(V))
  # mean(z) # etc.
  return(a)
}
