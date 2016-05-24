#' Supervised Component Generalized Linear Regression
#'
#' SCGLR implements a new Partial Least Squares regression approach in the multivariate generalized linear framework. The method
#' allows the joint modeling of random variables from different exponential family distributions, searching for 
#' common PLS-type components. \code{\link{scglr}} and \code{\link{scglrCrossVal}} are the two main functions.
#' The former constructs the components and performs the parameter estimation, while the
#' latter selects the approriate number of components by cross-validation. 
#' Dedicated plots, print, and summary functions are available.
#' The package contains also an ecological 
#' dataset dealing with the abundance of multiple tree genera given a large number of geo-referenced environmental
#' variables.
#' @references Bry X., Trottier C., Verron T. and Mortier F. (2013) Supervised Component Generalized Linear Regression using a PLS-extension of the Fisher scoring algorithm. \emph{Journal of Multivariate Analysis}, 119, 47-60.#' @docType package
#' @name scglr-package
#' @author Mortier F., Trottier C., Cornu G., Bry X.
#' @importFrom Matrix bdiag
#' @import Formula
#' @import expm
#' @importFrom graphics barplot
#' @import ggplot2
#' @import pROC
NULL

.onLoad <- function(libname,pkgname) {
}

.onAttach <- function(libname,pkgname) {
}
