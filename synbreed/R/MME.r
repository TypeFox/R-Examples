MME <- function(X,Z,GI,RI,y){
  SST <- t(y) %*% RI %*% y
  XX <- t(X) %*% RI %*% X
  XZ <- t(X) %*% RI %*% Z
  ZZ <- t(Z) %*% RI %*% Z
  Xy <- t(X) %*% RI %*% y
  Zy <- t(Z) %*% RI %*% y
  R1 <- cbind(XX,XZ)
  R2 <- cbind(t(XZ),(ZZ+GI))
  LHS <- rbind(R1,R2)
  RHS <- rbind(Xy,Zy)
  C <- ginv(LHS)
  bhat <- C %*% RHS
  SSR <- t(bhat) %*% RHS
  SSE <- SST - SSR
  N <- nrow(X)
  rX <- sum(diag(X%*%ginv(X)))
  sigma2e <- SSE/(N-rX)
  SEP <- sqrt(diag(C)*sigma2e)
  b <- bhat[1:ncol(X)]
  u <- bhat[-c(1:ncol(X))]
  return(list(b=b,u=u,LHS=LHS,RHS=RHS,C=C,SEP=SEP,SST=SST,SSR=SSR,residuals=y-X%*%b-Z%*%u))
  }
  
