lnL.M <-
function(D,
                  At,
                  at,
                  npar,
                  ret.alphas = FALSE)
{
  #Set up the D matrix
  Dm <- matrix(rep(0, npar^2), ncol = npar)
  Dm[upper.tri(Dm, diag = TRUE)] <- D
  #Make sure the diagonals of the D matrix are positive
  diag(Dm) <- abs(diag(Dm))
  #Estimate At + M and inv(At + M) for each year, where M = crossprod(D)
  AM <- mapply('+', At, list(crossprod(Dm)), SIMPLIFY = FALSE)
  AMi <- mapply(inv, AM, SIMPLIFY = FALSE)
  #Calculate AMi %*% at
  AMa <- mapply('%*%', AMi, at, SIMPLIFY = FALSE)
  #Calculate the sum of AMi
  AMis <- AMi[[1]]
  for (i in 2 : length(AMi))
    AMis <- AMis + AMi[[i]]
  #Calculate the sum of AMa %*% at
  AMas <- AMa[[1]]
  for (i in 2 : length(AMa))
    AMas <- AMas + AMa[[i]]
  #Calculate alphas given M
  aM <- inv(AMis) %*% AMas
  #If desired return alpha estimates and quit function, else continue
  #lnL calculation
  if(ret.alphas == TRUE)
  {
    #Add row and column names to the matrix showing which estimate
    #belongs to which alpha
    dimnames(aM) = list(names(at[[1]]), "")
    #Return alphas
    return(t(aM))
  }
  #Calculate the log likelihood function
  lnAM <- mapply(function(a) 2 * sum(log(diag(chol(a)))), AM, SIMPLIFY = FALSE)
  atmin <- mapply('-', at, list(aM), SIMPLIFY = FALSE)
  lnL <- sum(unlist(mapply(function(a, b, c){a + t(b) %*% c %*% b},
                           lnAM, atmin, AMi, SIMPLIFY = FALSE)))
  #Return
  lnL
}
