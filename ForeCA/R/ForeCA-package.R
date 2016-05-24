#' @title Implementation of Forecastable Component Analysis (ForeCA)
#' 
#' @description
#' Forecastable Component Analysis (ForeCA) is a novel dimension reduction
#' technique for multivariate time series \eqn{\mathbf{X}_t}.  
#' ForeCA finds a linar combination
#' \eqn{y_t = \mathbf{X}_t \mathbf{v}} that is easy to forecast. The measure of 
#' forecastability \eqn{\Omega(y_t)} (\code{\link{Omega}}) is based on the entropy 
#' of the spectral density \eqn{f_y(\lambda)} of \eqn{y_t}: higher entropy means 
#' less forecastable, lower entropy is more forecastable.
#' 
#' The main function \code{\link{foreca}} runs ForeCA on a 
#' multivariate time series  \eqn{\mathbf{X}_t}.
#'  
#' Please consult the \code{NEWS} file for a list of changes to previous
#' versions of this package.
#' 
#' @import MASS stats graphics reshape2 utils
#' @name ForeCA-package
#' @aliases ForeCA-package ForeCA
#' @docType package
#' @author Author and maintainer: Georg M. Goerg <im@@gmge.org>
#' @references 
#' Goerg, G. M. (2013). \dQuote{Forecastable Component Analysis}. 
#' Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
#' Available at \url{jmlr.org/proceedings/papers/v28/goerg13.html}.
#' @keywords package
#' @examples
#' XX <- ts(diff(log(EuStockMarkets)))
#' Omega(XX)
#' 
#' plot(log10(lynx))
#' Omega(log10(lynx))
#' 
#' \dontrun{
#' ff <- foreca(XX, n.comp = 4)
#' ff
#' plot(ff)
#' summary(ff)
#' }
#' 
NULL



