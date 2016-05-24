#' Generating selection gradients/vectors in random directions.
#' 
#' \code{randomBeta} generates unit length vectors (selection gradients) 
#' uniformly distributed in a k-dimensional hypersphere.
#' 
#' \code{randomBeta} exploits the spherical symmetry of a multidimensional 
#' Gaussian density function. Each element of each vector is randomly sampled 
#' from a univariate Gaussian distribution with zero mean and unit variance. The
#' vector is then divided by its norm to standardize it to unit length.
#' 
#' @param n Number of selection gradients/vectors.
#' @param k Number of dimensions.
#' @return \code{randomBeta} returns a matrix where the vectors are stacked column wise.
#' @author Geir H. Bolstad
#' 
#' @examples
#' ## Two vectors of dimension 3:
#' randomBeta(n = 2, k = 3)  
#' @export

randomBeta = function(n = 1, k = 2){
  X = matrix(rnorm(n*k), ncol = n)
  X = t(t(X)/sqrt(colSums(X^2)))
  rownames(X) = paste("dim", 1:k, sep="")
  X
}


