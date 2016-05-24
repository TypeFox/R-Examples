"Lcomoment.matrix" <-
function(DATAFRAME, k=1) {
  # DATAFRAME is data.frame of rectangular dimension
  # k is the kth order of L-comoments
  
  f <- length(DATAFRAME)        # how many fields or "random variables"
  M <- matrix(nrow=f, ncol=f)   # generate square matrix
  n <- length(DATAFRAME[,1])    # sample size 

  for(x1 in seq(1,f)) {         # BEGIN LOOP 1
    X1 <- DATAFRAME[,x1]        # extract array "1"
    for(x2 in seq(1,f)) {       # BEGIN LOOP 2
      X2 <- DATAFRAME[,x2]      # extract array "2"
      M[x1,x2] <- Lcomoment.Lk12(X1,X2,k) # compute the L-comoments
                                # for 1 and 2 and order k
      # If the moment order is first and the x1 and x2
      # indices are not equal, then off diagonals are NA
      if(k == 1 & x1 != x2)  M[x1,x2] = NA
    }                           # END LOOP 2
  }                             # END LOOP 1
  z <- list(type="Lcomoment.matrix", order=k, matrix=M)
  return(z)                     # return the matrix
}

