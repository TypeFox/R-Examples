"Lcomoment.coefficients" <-
function(Lk,L2) {
  # L2 is generally to denote a 2nd-order L-comoment matrix EXCEPT
  # when Lk is first order for L-CV like coefficients.

  # Following notation of Serfling and Xiao (2006)
  #  compute the L-comoment coefficients
  #  The univariate L-moment ratios are on the diagonal of Lk
  if(is.null(Lk$type) || Lk$type != "Lcomoment.matrix") {
    warning("First argument does not appear to be an L-comoment matrix")
    return()
  }
  if(is.null(L2$type) || L2$type != "Lcomoment.matrix") {
    warning("Second argument does not appear to be an L-comoment matrix")
    return()
  }
  if(Lk$order >= 2 && L2$order != 2) {
    warning("First L-comoment matrix is order 2 or greater, but second matrix is not of order 2")
    return()
  }
  if(Lk$order == 1) {
    if(L2$order != 2) { # In L-CV calculations L2/Lk, where Lk is 1st order,
                        # but in all others Lk/L2
      warning("First L-comoment matrix is order 1, but second matrix is not 2nd order.")
      return()
    }
    else {
      TMP <- Lk # swap Lk and L2 so that a single loop
      Lk  <- L2 # after this conditional can be used.
      L2  <- TMP
    }
  }
  LC      <- Lk$matrix        # to get the structure of Lk
  Lscales <- diag(L2$matrix)  # get univariate L-scale values
  n       <- length(Lscales)  # how many are there (how many columns)
  for(i in seq(1,n)) {        # loop through each column
    Lscale <- Lscales[i]      # extract single L-scale value
    LC[i,] <- Lk$matrix[i,]/Lscale   # divide the column by L-scale
                                     # to form coefficients
  }                           # end of loop
  z <- list(type="Lcomoment.coefficients", order = Lk$order, matrix = LC)
  return(z)                   # return the matrix
}

