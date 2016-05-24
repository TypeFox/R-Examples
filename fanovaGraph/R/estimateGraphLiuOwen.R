estimateGraphLiuOwen <- function(f.mat, d, q, q.arg, n.lo, confint, print.loop.index, ...) {
  ij.set <- combn(d, 2)
  p <- choose(d, 2)
  ZX01 <- matrix(runif(n.lo * 2*d), ncol = 2*d)
  ZX <- matrix(,n.lo,2*d)
  for (j in 1:(2*d)) ZX[, j] <- do.call(rep(q,2)[j], 
          c(list(p = ZX01[, j]),rep(q.arg,2)[[j]]))
  Z <- ZX[,1:d]
  X <- ZX[,-(1:d)]
  yZ <- f.mat(Z, ...)  
  y.i <- matrix(,n.lo,d)
  for (i in 1:d){
    if(print.loop.index) cat("index =",i ,"\n")
    Xi <- Z
    Xi[,i] <- X[,i]
    y.i[,i] <- f.mat(Xi, ...)
  }
  y.ij <- matrix(,n.lo,p)
  for (r in 1:p){
    ij <- ij.set[,r]
    if(print.loop.index) cat("index =",ij,"\n")
    Xij <- Z
    Xij[,ij] <- X[,ij]
    y.ij[,r] <- f.mat(Xij, ...)
  }
  totalInt <- numeric(p)
  if (confint) {
    varInt <- numeric(p)
  }
  for (index in 1:p)
  {
    i <- ij.set[1,index]; 
    j <- ij.set[2,index]
    deltasq<-(y.ij[,index]-y.i[,i]-y.i[,j]+yZ)^2
    totalInt[index] <- 1/(4*n.lo) * sum(deltasq)
    if (confint) {
      varInt[index] <- var(deltasq) / 16
    }
  }    
  inter <- paste("X", ij.set[1, ], "*", "X", ij.set[2,], sep = "")
  res <- as.matrix(totalInt)
  if (confint){ 
    Std.Error <- sqrt(varInt)/sqrt(n.lo)
    lower <- totalInt - qnorm(0.975)*Std.Error
    upper <- totalInt + qnorm(0.975)*Std.Error
    res <- cbind(totalInt, Std.Error,lower,upper)
  }
  rownames(res) <- inter
  return(res)
}