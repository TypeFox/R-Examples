#####################
#' Robust covariance matrix estimation
#' @description Obtains a robust estimate of the covariance matrix of a sample of multivariate data, 
#'              using Campbell's (1980) method as described on p231-235 of Krzanowski (1988). 
#' @param sY A matrix, where each column is a replicate observation on a multivariate r.v.
#' @param alpha tuning parameter, see details.
#' @param beta tuning parameter, see details.
#' @details Campbell (1980) suggests an estimator of the covariance matrix which downweights observations 
#'                          at more than some Mahalanobis distance \code{d.0} from the mean.
#'                          \code{d.0} is \code{sqrt(nrow(sY))+alpha/sqrt(2)}. Weights are one for observations 
#'                          with Mahalanobis distance, \code{d}, less than \code{d.0}. Otherwise weights are 
#'                          \code{d.0*exp(-.5*(d-d.0)^2/beta)/d}. The defaults are as recommended by Campbell.
#'                          This routine also uses pre-conditioning to ensure good scaling and stable 
#'                          numerical calculations.
#' @return A list where:
#'         \itemize{
#'         \item{\code{E}}{a square root of the inverse covariance matrix. i.e. the inverse cov 
#'                         matrix is \code{t(E)\%*\%E};}
#'         \item{\code{half.ldet.V}}{Half the log of the determinant of the covariance matrix;}
#'         \item{\code{mY}}{The estimated mean;} 
#'         \item{\code{sd}}{The estimated standard deviations of each variable.}
#'          }
#' @references Krzanowski, W.J. (1988) Principles of Multivariate Analysis. Oxford.
#'             Campbell, N.A. (1980) Robust procedures in multivariate analysis I: robust covariance estimation. JRSSC 29, 231-237. 
#' @author Simon N. Wood, maintained by Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' p <- 5;n <- 100
#' Y <- matrix(runif(p*n),p,n)
#' robCov(Y)
#' @export

robCov <- function(sY, alpha=2, beta=1.25) {

  .robCov(sY = sY, 
          alpha = alpha, 
          beta = beta, 
          alpha2 = NULL, 
          beta2 = NULL, 
          tolVar = 0.0)
  
}
