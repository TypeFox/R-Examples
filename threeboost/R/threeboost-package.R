#' @title
#' Thresholded boosting based on estimating equations
#'
#' @description
#' The \code{threeboost} package implements a the EEBoost algorithm described in \emph{Wolfson (2011, JASA)}.
#' EEBoost is a general-purpose method for variable selection which can be applied whenever inference would be based on an estimating equation.
#' Thresholded EEBoost (function \code{\link{threeboost}}) is a generalization of EEBoost which allows multiple variables to enter the model at each boosting step. 
#' EEBoost (function \code{\link{eeboost}}) is a special case of thresholded boosting with the threshold set to 1.
#' 
#' The package currently provides a "pre-packaged" function, \code{\link{geeboost}}, which carries out
#' variable selection for correlated outcome data based on the Generalized Estimating Equations. However, 
#' the \code{\link{threeboost}} (and \code{\link{eeboost}}) functions can also accommodate user-provided estimating functions. 
#'
#' @docType package
#' @name threeboost-package
NULL