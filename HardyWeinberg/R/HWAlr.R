HWAlr <- function(X,zeroadj=0.5,denominator=2) {
  # compute the additive logratio transformation for a vector or for each row of a matrix
  X[X==0] <- zeroadj
  numerator <- 1:3
  numerator <- setdiff(numerator,denominator)
  if(is.matrix(X)) {
     c1 <- log(X[,numerator[1]]/X[,denominator])
     c2 <- log(X[,numerator[2]]/X[,denominator])
     Y <- cbind(c1,c2)
   }
  if(is.vector(X)) {
    c1 <- log(X[numerator[1]]/X[denominator])
    c2 <- log(X[numerator[2]]/X[denominator])
    Y <- c(c1,c2)
  }
  return(Y)
}
