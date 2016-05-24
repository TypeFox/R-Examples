##############################################################
#' Fast simulations of multivariate normal random variables
#'
#' @param n number of random vectors to be simulated.
#' @param mu vector of length d, representing the central position.
#' @param sigma covariance matrix (d x d). Alternatively is can be the cholesky decomposition
#'              of the covariance. In that case \code{isChol} should be set to \code{TRUE}.
#' @param ncores Number of cores used. The parallelization will take place only if OpenMP is supported.
#' @param isChol boolean set to true is \code{sigma} is the cholesky decomposition of the covariance matrix.
#' @return A vector of length n where the i-the entry contains the square mahalanobis distance i-th random vector.
#' @details Notice that this function does not use one of the Random Number Generators (RNGs) provided by R, but one 
#'          of the parallel cryptographic RNGs described in (Salmon et al., 2011). It is important to point out that this
#'          RNG can safely be used in parallel, without risk of collisions between parallel sequence of random numbers.
#'          The initialization of the RNG depends on R's seed, hence the \code{set.seed()} function can be used to 
#'          obtain reproducible results. Notice though that changing \code{ncores} causes most of the generated numbers
#'          to be different even if R's seed is the same (see example below). NB: at the moment the RNG does not work
#'          properly on Solaris OS when \code{ncores>1}. Hence, \code{rmvn()} checks if the OS is Solaris and, if this the case, 
#'          it imposes \code{ncores==1}. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>, C++ RNG engine by Thijs van den Berg <http://sitmo.com/>.
#' @references  John K. Salmon, Mark A. Moraes, Ron O. Dror, and David E. Shaw (2011). Parallel Random Numbers: As Easy as 1, 2, 3.
#'              D. E. Shaw Research, New York, NY 10036, USA.
#' @examples
#' d <- 5
#' mu <- 1:d
#' 
#' # Creating covariance matrix
#' tmp <- matrix(rnorm(d^2), d, d)
#' mcov <- tcrossprod(tmp, tmp)
#' 
#' set.seed(414)
#' rmvn(4, 1:d, mcov)
#' 
#' set.seed(414)
#' rmvn(4, 1:d, mcov)
#' 
#' set.seed(414)  
#' rmvn(4, 1:d, mcov, ncores = 2) # r.v. generated on the second core are different
#' 
#' @export rmvn

rmvn <- function(n, mu, sigma, ncores = 1, isChol = FALSE)
{
  
  if( !is.matrix(sigma) ) sigma <- as.matrix( sigma )
  
  if( ncores > 1 && grepl('SunOS', Sys.info()['sysname']) ){
    
    message("rmvn() cannot be used on multiple cores under Solaris. I am resetting \"ncores\" to 1.")
    
    ncores <- 1
    
  }

  .Call( "rmvnCpp", 
         n_ = n, 
         mu_ = mu, 
         sigma_ = sigma, 
         ncores_ = ncores,
         isChol_ = isChol, 
         PACKAGE = "mvnfast" )
  
}
