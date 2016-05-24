NULL

#' 
#' \code{GPCAiteration} S3 class returned by \code{\link{GPCA_iteration}}
#' 
#' 	\describe{
#' 	\item{\code{x_prev}}{ Previous set of random variable, \code{x_prev} input variable of  \code{\link{GPCA_iteration}} }
#' 
#' \item{\code{x_gauss_prev}}{ Marginal Gaussianization of \code{x_prev} obtained through \code{\link{normalizeGaussian_severalstations}}}
#' 
#' \item{\code{B_prev}}{ rotation matrix (i. e. eigenvector matrix of the covariance matrix of  \code{x_gauss_prev})}
#' 
#' \item{\code{x_next} }{results obtained by multiplying \code{B_prev} by  \code{x_gauss_prev} (see equation 1 of  the reference in \code{\link{GPCA_iteration}}) }
#'  }
#' 
#' 
#' @title GPCAiteration-class 
#' 
#' @note Formal definition with \code{\link{setOldClass}} for the S3 class \code{GPCAiteration}
#' 
#' 

#' 
#' 
#' @author Emanuele Cordano
#' 
#' @docType class
#' @aliases GPCAiteration
#' @name GPCAiteration-class
#' @rdname GPCAiteration-class
#' 
#' @keywords classes
#' @exportClass GPCAiteration
#'
#' @examples showClass("GPCAiteration")
#' 
#'  
#' 
setOldClass("GPCAiteration")
