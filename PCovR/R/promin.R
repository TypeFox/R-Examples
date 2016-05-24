promin <-
function(F1,nrs=20){
  r1 <- wvarim(F1,nrs)
  Tv <- r1$Th
  F <- F1 %*% Tv
  
  #normalize
  A <- diag(diag(F %*% t(F)))
  c <- eigen(A)$vectors%*%diag(1/sqrt(eigen(A)$values))%*%t(eigen(A)$vectors)
  c <- c %*% F
  
  #majorization of a vector for promin
  PHI <- c^2
  n <- nrow(PHI)
  j <- ncol(PHI)
  m <- colMeans(PHI)
  s <-apply(PHI,2,sd)*sqrt((n-1)/n)
  c <- m + s/4
  cvec <- t(matrix(c,j,n))
  PHIref <- PHI-cvec
  B <- PHI*(PHIref>=0)
  W <- B>0.0001
  
  #semi-specified oblique rotation
  r2 <- tarrotob(F1, W)
  T <- r2$T
  U <- solve(t(T))
  loadings <- F1 %*% T
  results <- list(Th=T,loadings=loadings,U=U)
  return(results)
}
