#' Compute the mean square error between \eqn{X_t} and \eqn{Y_t}, as \deqn{\frac{1}{n} \sum_{k=1}^n \|X_k-Y_k\|^2,} 
#' where \eqn{\|\cdot\|} denotes a Euclidian norm and \eqn{n} is the number of observations.
#'
#' @title Compute a mean square error between X and Y
#' @param X first matrix to compare
#' @param Y second matrix to compare
#' @return Estimated mean square error
#' @export 
MSE = function(X,Y){ 
  if (is.vector(X))
    X = matrix(X)
  if (is.vector(Y) && Y!= 0)
    Y = matrix(Y)

  if (!is.matrix(X) || (!is.matrix(Y) && Y != 0))
    stop("X and Y must be matrices, or Y=0")
  if (any(dim(Y) != dim(X)) && Y != 0)
    stop("Dimentions of X and Y must be equal")
  
  S = 0
  for (i in 1:dim(X)[1]){
    S = S + abs((Conj((Y-X)[i,])) %*%(Y-X)[i,])
  }
  S/dim(X)[1]
}