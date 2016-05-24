vifx <- function(X){
  Xs <- scale(X)
  xx <- t(Xs)%*%Xs/(nrow(X)-1) #put in correlation form
  diag(solve(xx)) #VIF are the diagnonal elements of the inverse
}

