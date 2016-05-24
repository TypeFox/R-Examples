#' Function \code{tailIntEmp}
#' 
#' Integral of the bivariate empirical stable tail dependence function over the unit square.
#' 
#' @param ranks A \code{n} x 2 matrix, where each column is a permutation of the integers \code{1:n}, representing the ranks computed from a sample of size \code{n}.
#' @param n The sample size. If not specified, it is computed as the number of rows of the matrix \code{ranks}.
#' @param k The threshold parameter in the definition of the empirical stable tail dependence function. An integer between 1 and \code{n-1}.
#' @return A scalar.
#' @details This is an analytic implementation of the integral of the stable tail dependence function,
#' which is much faster than numerical integration.
#' See Einmahl et al. (2014) for a definition of the empirical stable tail dependence function. 
#' @references Einmahl, J.H.J., Kiriliouk, A., Krajina, A. and Segers, J. (2014), "An M-estimator of spatial tail dependence". See \url{http://arxiv.org/abs/1403.1975}. 
#' @export
#' @seealso \code{\link{Mestimator}}, \code{\link{tailInt}}
#' @examples
#' n <- 20
#' (ranks <- cbind(sample(1:n), sample(1:n)))
#' tailIntEmp(ranks, k = 5)
tailIntEmp <- function(ranks, n = nrow(ranks), k) {
  res<-.C("stdf",as.integer(ranks[, 1]),as.integer(ranks[, 2]),as.integer(n),
          as.integer(k),result=double(1),PACKAGE="spatialTailDep")$result
  return(res)
}