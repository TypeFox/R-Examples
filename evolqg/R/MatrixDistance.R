#' Matrix distance
#' 
#' Calculates Distances between covariance matrices.
#' 
#' @param cov.x Single covariance matrix or list of covariance matrices.
#' If single matrix is suplied, it is compared to cov.y.
#' If list is suplied and no cov.y is suplied, all matrices
#' are compared.
#' If cov.y is suplied, all matrices in list are compared to it.
#' @param cov.y First argument is compared to cov.y.
#' Optional if cov.x is a list.
#' @param distance distance function for use in calculation. Currently supports "Riemann" and "Overlap".
#' @param ... aditional arguments passed to other methods
#' @param parallel if TRUE and a list is passed, computations are done in parallel. Some foreach backend must be registered, like doParallel or doMC.
#' @return
#' If cov.x and cov.y are passed, returns distance between them.
#'
#' If is a list cov.x and cov.y are passed, same as above, but for all matrices in cov.x.
#'
#' If only a list is passed to cov.x, a matrix of Distances is returned
#' @export
#' @rdname MatrixDistance
#' @author Diogo Melo
#' @seealso \code{\link{RiemannDist}},\code{\link{OverlapDist}}
#' @examples
#' c1 <- RandomMatrix(10)
#' c2 <- RandomMatrix(10)
#' c3 <- RandomMatrix(10)
#' MatrixDistance(c1, c2, "OverlapDist")
#' MatrixDistance(c1, c2, "RiemannDist")
#'
#' MatrixDistance(list(c1, c2, c3), distance = "OverlapDist")
#'
#'
#' c4 <- RandomMatrix(10)
#' MatrixDistance(list(c1, c2, c3), c4)
#' @keywords matrixcomparison
#' @keywords matrixDistance
MatrixDistance <- function(cov.x, cov.y, distance, ...) 
  UseMethod("MatrixDistance")

#' @rdname MatrixDistance
#' @method MatrixDistance default
#' @export
MatrixDistance.default <- function (cov.x, cov.y, distance = c('OverlapDist', 'RiemannDist'), ...) {
  distance <- match.fun(match.arg(distance))
  output <- distance(cov.x, cov.y)
  return(output)
}

#' @rdname MatrixDistance
#' @method MatrixDistance list
#' @export
MatrixDistance.list <- function (cov.x, cov.y = NULL, distance = c('OverlapDist', 'RiemannDist'), ...,  parallel = FALSE)
{
  distance <- match.fun(match.arg(distance))
  if (is.null (cov.y)) {
    output <- ComparisonMap(cov.x,
                         function(x, y) return(c(distance(x, y,...), NA)),
                         parallel = parallel)
    output <- output[[1]]
  } else{
    output <- SingleComparisonMap(cov.x, cov.y,
                               function(x, y) return(c(distance(x, y,...), NA)),
                               parallel = parallel)
    output <- output[-length(output)]
  }
  return(output)
}


#' Matrix Riemann distance
#'
#' Return distance between two covariance matrices
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @return Riemann distance between cov.x and cov.y
#' @author Edgar Zanella
#' @export
#' @references Mitteroecker, P., & Bookstein, F. (2009). The ontogenetic trajectory of the phenotypic covariance matrix, with examples from craniofacial shape in rats and humans. Evolution, 63(3), 727-737. doi:10.1111/j.1558-5646.2008.00587.x
#' @keywords matrixDistance
#' @keywords matrixcomparison

RiemannDist <- function(cov.x, cov.y) {
  return (sqrt(sum(log(eigen(solve(cov.x, cov.y))$values)^2)))
}

#'Distribution overlap distance
#'
#'Calculates the overlap between two normal distributions, 
#'defined as the probability that a draw from one distribution 
#'comes from the other
#'
#' @param cov.x covariance or correlation matrix
#' @param cov.y covariance or correlation matrix
#' @param iterations number of drows
#' @return Overlap distance between cov.x and cov.y
#' @references Ovaskainen, O. (2008). A Bayesian framework for comparative quantitative genetics. Proceedings of the Royal Society B, 669-678. doi:10.1098/rspb.2007.0949
#' @export
#' @importFrom mvtnorm rmvnorm dmvnorm
OverlapDist <- function(cov.x, cov.y, iterations = 10000){
  x_samples <- rmvnorm(iterations, sigma = cov.x)
  y_density <- dmvnorm(x_samples, sigma = cov.y)
  x_density <- dmvnorm(x_samples, sigma = cov.x)
  overlap <- mean(y_density/(x_density + y_density))  
  return(sqrt(1-overlap*2))
}
