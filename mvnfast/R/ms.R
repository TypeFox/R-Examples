#####################
### Mean-shift algorithm
#####################
#####
#' Mean-shift mode seeking algorithm
#' @description Given a sample from a d-dimensional distribution, an initialization point and a bandwidth
#'              the algorithm finds the nearest mode of the corresponding Gaussian kernel density.
#' @param X n by d matrix containing the data.
#' @param init d-dimensional vector containing the initial point for the optimization. By default
#'             it is equal to \code{colMeans(X)}.
#' @param H Positive definite bandwidth matrix representing the covariance of each component of the Gaussian kernel density.
#' @param tol Tolerance used to assess the convergence of the algorithm, which is stopped if the absolute values
#'            of increments along all the dimensions are smaller then tol at any iteration. Default value is 1e-6.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param store If \code{FALSE} only the latest iteration is returned, if \code{TRUE} the function will return a matrix where
#'             the i-th row is the position of the algorithms at the i-th iteration.
#' @return A list where \code{estim} is a d-dimensional vector containing the last position of the algorithm, while \code{traj} 
#'         is a matrix with d-colums representing the trajectory of the algorithm along each dimension. If \code{store == FALSE} the whole trajectory
#'         is not stored and \code{traj = NULL}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' set.seed(434)
#' 
#' # Simulating multivariate normal data
#' N <- 1000
#' mu <- c(1, 2)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' X <- rmvn(N, mu = mu, sigma = sigma)
#' 
#' # Plotting the true density function
#' steps <- 100
#' range1 <- seq(min(X[ , 1]), max(X[ , 1]), length.out = steps)
#' range2 <- seq(min(X[ , 2]), max(X[ , 2]), length.out = steps)
#' grid <- expand.grid(range1, range2)
#' vals <- dmvn(as.matrix(grid), mu, sigma)
#' 
#' contour(z = matrix(vals, steps, steps),  x = range1, y = range2, xlab = "X1", ylab = "X2")
#' points(X[ , 1], X[ , 2], pch = '.')
#'  
#' # Estimating the mode from "nrep" starting points
#' nrep <- 10
#' index <- sample(1:N, nrep)
#' for(ii in 1:nrep) {
#'   start <- X[index[ii], ]
#'   out <- ms(X, init = start, H = 0.1 * sigma, store = TRUE)
#'   lines(out$traj[ , 1], out$traj[ , 2], col = 2, lwd = 2) 
#'   points(out$final[1], out$final[2], col = 4, pch = 3, lwd = 3) # Estimated mode (blue)
#'   points(start[1], start[2], col = 2, pch = 3, lwd = 3)         # ii-th starting value 
#' }
#' @export ms

ms <- function(X, init, H, tol = 1e-6, ncores = 1, store = FALSE)
{
  if( is.matrix(X) == FALSE ) X <- matrix(X, length(X), 1)
  
  if( !is.matrix(H) ) H <- diag(H, ncol(X))

  if( !store ) trajectory <- NULL
  
  cholDec <- chol( H )
  
  tmp <- .Call( "msCpp", 
                init_ = init,
                X_ = X, 
                cholDec_ = cholDec, 
                ncores_ = ncores,
                tol_ = tol, 
                store_ = store,
                PACKAGE = "mvnfast" )
  
  list("final" = drop(tmp$final), "traj" = do.call("rbind", tmp$traj))
}
