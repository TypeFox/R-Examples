HWClr <- function(X,zeroadj=0.5) {
  # compute the centred log-ratio transformation for each row of X.
  X[X==0] <- zeroadj
  if(is.matrix(X)) {
    gm <- apply(X,1,gmeanrow)
    Y <- log(solve(diag(gm))%*%X)
  }
  if(is.vector(X)) {
    gm <- gmeanrow(X)
    Y <- log((1/gm)*X)
  }
  return(Y)
}
