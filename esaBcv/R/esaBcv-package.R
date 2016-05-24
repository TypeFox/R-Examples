#' esaBcv
#'
#' The esaBcv package provides functions to estimate the latent factors of a given 
#' matrix, no matter it is high-dimensional or not. It tries to first estimate the 
#' number of factors using Bi-cross-validation and then estimate the latent factor 
#' matrix and the noise variances using an Early-stopping-alternation method.
#' The method is proposed by Art B. Owen and Jingshu Wang (2015). 
#' @name esaBcv_package
#' @docType package
#' @author Maintainer: Jingshu Wang <wangjingshususan@@gmail.com>
#' @seealso Owen and Wang (2015) Bi-cross-validation for factor analysis,
#' \url{http://arxiv.org/abs/1503.03515}
#' @examples
#'
#' \dontrun{
#' data(simdat)
#' result <- EsaBcv(simdat$Y)
#' plot(result)
#' 
#' }
#'
NULL
