nfs <-
function(At,
                at,
                npar,
                nyear)
{
  #Calculate inv(At)
  Ati <- lapply(At, inv)
  #Calculate the covariance matrix under the assumption of no
  #fluctuation selection
  ret <- list(Anf = Ati[[1]])
  for (i in 2 : length(Ati))
    ret$Anf <- ret$Anf + Ati[[i]]
  ret$Anf <- inv(ret$Anf)
  #Add row and column names to Anf
  dimnames(ret$Anf) <- dimnames(At[[1]])
  #Calculate Ati multiplied with at
  Atia <- mapply('%*%', Ati, at, SIMPLIFY = FALSE)
  #Sum Atia over years (t)
  Aia <- Atia[[1]]
  for ( i in 2 : length(Atia))
    Aia <- Aia + Atia[[i]]
  #Calculate the alpha vector under the assumption of no
  #fluctuation selection
  ret$anf <- ret$Anf %*% Aia
  #Add row and column names to the matrix showing which estimate belongs
  #to which alpha
  dimnames(ret$anf) <- list(names(at[[1]]), NULL)
  ret$anf <- t(ret$anf)
  #Output
  ret
}
