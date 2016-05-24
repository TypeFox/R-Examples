#-----------------------------------------------------------------------------------------------------------------#
# Function SVMLinearMonotone
#-----------------------------------------------------------------------------------------------------------------#

#Solve the linear SVM problem with monotonicity constraints
#' @import quadprog
SVMLinearMonotone <- function(X,Y){
  ndim <- dim(X)[2]
  dm   <- matrix(0, ndim + 1)
  Am   <- rbind(cbind( matrix(rep(Y, ndim), ncol = ndim)*X, Y), cbind(diag(ndim), rep(0, ndim)))
  b_0m <- rbind( matrix(1, nrow = dim(X)[1]), matrix(0, nrow = ndim)) 
  Dm   <- diag(ndim + 1)
  Sm   <- solve.QP(Dmat = Dm, dvec = dm, Amat = t(Am), bvec = b_0m)
  return(Sm$solution)
}


