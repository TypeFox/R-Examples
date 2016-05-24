totalIndex <- function(f.mat, d, q=NULL, q.arg=NULL, n.mc, ...)
{
  if (is.null(q)) {     
    q <- rep("qunif", d)
  } else if (length(q) == 1) {
    q <- rep(q, d)
  }
  if (is.null(q.arg)) {
    q.arg <- rep(list(list()), d)
  } else if (FALSE %in% sapply(q.arg, is.list)) {
    q.arg <- rep(list(q.arg), d)
  }
  X <- matrix(runif(n.mc * d), ncol = d)
  for (j in 1:d) X[, j] <- do.call(q[j], c(list(p = X[,j]), q.arg[[j]]))
  Z <- matrix(runif(n.mc * d), ncol = d)
  for (j in 1:d) Z[, j] <- do.call(q[j], c(list(p = Z[,j]), q.arg[[j]]))
  yZ <- f.mat(Z,...)
  
  DT <- numeric(d)
  names(DT) <- 1:d
  for (i in 1:d){
    ZXi <- Z
    ZXi[,i] <- X[,i]
    DT[i] <- 1/(n.mc*2)*sum( (f.mat(ZXi,...) - yZ)^2)
  }
  return(DT)
}
