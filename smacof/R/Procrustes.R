Procrustes <- function(X, Y){
  n <- dim(X)[1]
  E <- diag(1, nrow = n)
  eins <- rep(1, n)
  k <- 1/n
  Z <- E-k*eins %*% t(eins)
  C <- t(X)%*%Z%*%Y
  s <- svd(C) 
  f <- diag(s$d)
  P <- s$u
  Q <- s$v
  T <- Q%*%t(P)
  streck <- sum(diag((t(X) %*% Z %*% Y %*% T)))/sum(diag((t(Y) %*% Z %*% Y)))
  trans <- as.vector(k * t(X-streck*Y %*% T) %*% eins)
  Yhut <- streck*Y%*%T + eins%*%t(trans)
  
  colnames(Yhut) <- rownames(T) <- colnames(T) <- names(trans) <- colnames(Y)
  
  ## distances and congruence coefficient
  dX <- dist(X)
  dY <- dist(Y)
  dYhat <- dist(Yhut)
  cong <- sum(dX*dY)/(sqrt(sum(dX^2))*sqrt(sum(dY^2)))
  
  res <- list(X = X, Y = Y, Yhat = Yhut, translation = trans, dilation = streck, rotation = T, 
              confdistX = dX, confdistY = dY, confdistYhat = dYhat, congcoef = cong, call = match.call())
  class(res) <- "procr"
  return(res)
}