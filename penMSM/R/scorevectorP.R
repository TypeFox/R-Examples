scorevectorP <- function(mu, X, event){
  return(as.numeric(crossprod(X, event - mu)))
  ## as.numeric(t(X) %*% (event - mu))
}

