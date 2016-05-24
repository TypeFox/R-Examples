sF <- function(mu, X, event){
  return(list(s = c(crossprod(X, event - mu)),
              F = crossprod(X, diag(mu)) %*% X))
}