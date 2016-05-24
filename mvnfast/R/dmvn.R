#########
#### Fast computation of multivariate normal pdf
#########
#' Fast computation of the multivariate normal density.
#'
#' @param X matrix n by d where each row is a d dimensional random vector. Alternatively \code{X} can be a d-dimensional vector.
#' @param mu vector of length d, representing the mean the distribution.
#' @param sigma covariance matrix (d x d). Alternatively is can be the cholesky decomposition
#'              of the covariance. In that case isChol should be set to TRUE.
#' @param log boolean set to true the logarithm of the pdf is required.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param isChol boolean set to true is \code{sigma} is the cholesky decomposition of the covariance matrix.
#' @return a vector of length n where the i-the entry contains the pdf of the i-th random vector.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> 
#' @examples
#' N <- 100
#' d <- 5
#' mu <- 1:d
#' X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
#' tmp <- matrix(rnorm(d^2), d, d)
#' mcov <- tcrossprod(tmp, tmp)
#' myChol <- chol(mcov)
#' 
#' head(dmvn(X, mu, mcov), 10)
#' head(dmvn(X, mu, myChol, isChol = TRUE), 10)
#' 
#' \dontrun{
#' # Performance comparison
#' library(mvtnorm)
#' library(microbenchmark)
#' 
#' a <- cbind(
#'       dmvn(X, mu, mcov),
#'       dmvn(X, mu, myChol, isChol = TRUE),
#'       dmvnorm(X, mu, mcov))
#'       
#' # Check if we get the same output as dmvnorm()
#' a[ , 1] / a[, 3]
#' a[ , 2] / a[, 3]
#' 
#' microbenchmark(dmvn(X, mu, myChol, isChol = TRUE), 
#'                dmvn(X, mu, mcov), 
#'                dmvnorm(X, mu, mcov))
#' }
#' @export dmvn

dmvn <- function(X, mu, sigma, log = FALSE, ncores = 1, isChol = FALSE){
  
  if( !is.matrix(X) ) X <- matrix(X, 1, length(X))
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
    
  .Call( "dmvnCpp", 
         X_ = X, 
         mu_ = mu, 
         sigma_ = sigma, 
         log_ = log, 
         ncores_ = ncores,
         isChol_ = isChol, 
         PACKAGE = "mvnfast" )
}
