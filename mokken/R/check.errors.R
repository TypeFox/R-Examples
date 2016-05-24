"check.errors" <- function(X){

  X <- check.data(X)
  maxx <- max(X)
  minx <- min(X)
  N <- nrow(X)
  J <- ncol(X)
  Y <- matrix(t(X),1,N*J)
  Z <- matrix(rep(Y,maxx),maxx,N*J,T)
  Z <- ifelse(Z < row(Z),0,1)
  Z <- matrix(as.vector(Z),N,maxx*J,T)

  # COMPUTE WEIGHTS

  if (maxx == 1) tmp.1 <- matrix(apply(X,2,tabulate, maxx), nrow=1) else tmp.1 <- apply(X,2,tabulate, maxx)
  tmp.2 <- apply(tmp.1,2,function(x) rev(cumsum(rev(x))))+runif(J*maxx,0,1e-3)
  # runif is added to avoid equal ranks
  tmp.3 <- matrix(rank(-tmp.2),1,maxx*J)
  # tmp.3 is a vector with the order of the ISRFs
  Z <- Z[,order(tmp.3)]
  w <- apply(Z,1,function(x){sum(x*cumsum(abs(x-1)))})
  return(w)
}
