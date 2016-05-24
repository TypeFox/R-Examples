#' The walkr package.
#' 
#' \tabular{ll}{ Package \tab walkr \cr Type: \tab Package\cr 
#' Version: \tab 0.3.1\cr Date: \tab 2015-07-14\cr License: \tab GPL-3\cr }
#'
#' The \code{walkr} package samples points using MCMC random walks from the 
#' intersection of the \eqn{N}-Simplex with \eqn{M} hyperplanes. Mathematically
#' speaking, the sampling space is all vectors \eqn{x} that satisfy
#' \eqn{Ax = b}, \eqn{\sum x  = 1}, and \eqn{x_i \ge 0}. The sampling
#' algorithms implemented are hit-and-run and Dikin Walk. \code{walkr} also 
#' provides tools to examine and visualize the convergence properties of the 
#' random walks.
#' 
#' The main function of the package is \code{walkr}. The user specifies \eqn{A} and
#' \eqn{b} in \eqn{Ax = b}, and the \code{walkr} function samples points 
#' from the complete solution to \eqn{Ax=b} intersected with the \eqn{N}-simplex. 
#' The user can choose either \code{"dikin"} or \code{"hit-and-run"} as the 
#' sampling method, and the function also provides other MCMC parameters 
#' such as thinning and burning. 
#' 
#' @name walkr
#' @docType package
#' 
#' @importFrom hitandrun har
#' @importFrom limSolve lsei
#' @importFrom MASS Null
#' @importFrom stats rnorm runif
#' @importFrom Rcpp evalCpp
#' 
#' @useDynLib walkr
#' @exportPattern "^[[:alpha:]]+"
NULL

