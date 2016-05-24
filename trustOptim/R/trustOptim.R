#' @title Trust-region optimization
#' 
#' @description Nonlinear optimizers using trust regions, with methods
#'   optimized for sparse Hessians.
#' 
#' @details
#' Trust region algorithm for nonlinear optimization. In addition to being more stable and robust than optim, this package includes methods that are scalable and efficient (in terms of both speed and memory usage) when the Hessian is sparse.
#' 
#' @references
#' Braun, Michael.  2014.  trustOptim:  An R Package for Trust Region
#' Optimization with Sparse Hessians. Journal of Statistical Software 60(4),
#' 1-16. www.jstatsoft.org/v60/i04/.
#' 
#' Nocedal, Jorge, and Stephen J Wright. 2006. Numerical Optimization.
#' Second edition. Springer.
#' 
#' Steihaug, Trond. 1983. The Conjugate Gradient Method and Trust
#' Regions in Large Scale Optimization. SIAM Journal on Numerical
#' Analysis 20(3), 626-637.
#'
#' @docType package
#' @name trustOptim
#' @import Matrix Rcpp
#' @useDynLib trustOptim
NULL 
