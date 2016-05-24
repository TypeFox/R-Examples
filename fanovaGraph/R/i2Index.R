i2Index <- function(f.mat, d, q=NULL, q.arg=NULL, n.i2, ...)
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
  ij.set <- combn(d, 2)
  p <- choose(d, 2)
  
  ZX01 <- matrix(runif(n.i2 * 2*d), ncol = 2*d)
  ZX <- matrix(,n.i2,2*d)
  for (j in 1:(2*d)) ZX[, j] <- do.call(rep(q,2)[j], c(list(p = ZX01[, j]),rep(q.arg,2)[[j]]))
  Z <- data.frame(ZX[,1:d])
  X <- data.frame(ZX[,-(1:d)])
  i2 <- as.matrix(sobol(model = f.mat, X1 = X, X2 = Z, order = 2, ...)$D[(d+1):(d+p),])
  rownames(i2) <- paste("X", ij.set[1, ], "*", "X", ij.set[2,], sep = "")  
  colnames(i2) <- paste("i2Index")
  return(i2)
}
