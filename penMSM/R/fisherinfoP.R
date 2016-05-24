fisherinfoP <- function(mu, X){
  return(crossprod(X, diag(mu)) %*% X)
  ## return(t(X) %*% Dbeta %*% X)
}