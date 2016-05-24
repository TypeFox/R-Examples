## Create a covariance matrix
S <- cbind(c(2,1),c(1,2))
## compute Cholesky factor
R <- chol(S)
## compute determinant
log(det(R))
## compare with sum of the logarithm of diagonal elements
sumLogDiag(R)
##or using sumLog (usefull e.g. for the Matrix-class)
sumLog(diag(R))

\dontshow{
  if( abs(log(det(R)) - sumLogDiag(R)) > 1e-10 ){
    stop("sumLogDiag: Results not equal")
  }
  if( abs(sumLog(diag(R)) - sumLogDiag(R)) > 1e-10 ){
    stop("sumLog: Results not equal")
  }
}

