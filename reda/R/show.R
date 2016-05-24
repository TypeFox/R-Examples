################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
#' @include class.R 
NULL


#' Show an object.
#' 
#' An S4 class generic function that displays certain object.
#' 
#' \itemize{
#'   \item For \code{\link{rateReg-class}} object, 
#'       it prints out brief summary of the fitted model.
#'   \item For \code{\link{summaryRateReg-class}} object, 
#'       it prints out summary of the fitted model.
#'   \item For \code{\link{sampleMcf-class}} object,
#'       it prints out the function call, formula and
#'       the sample MCF data frame.
#'   \item For \code{\link{rateRegMcf-class}} object,
#'       it prints formula, new data, confidence level,
#'       and the estimated MCF data frame.
#' }
#' 
#' @param object An object used to dispatch a method.
#' @name show-method
#' @seealso
#' \code{\link{rateReg}} for model fitting;
#' \code{\link{summary,rateReg-method}} for summary of a fitted model;
#' \code{\link{mcf}} for estimation of MCF.
#' @importFrom methods show
NULL


#' @rdname show-method
#' @aliases show,rateReg-method
#' @export
setMethod(f = "show", signature = "rateReg",
          definition = function(object) {
              beta <- object@estimates$beta[, 1]
              names(beta) <- row.names(object@estimates$beta)
              theta <- object@estimates$theta[, 1]
              names(theta) <- NULL
              alpha <- object@estimates$alpha[, 1]
              names(alpha) <- row.names(object@estimates$alpha)
              cat("Call: \n")
              print(object@call)
              cat("\nCoefficients of covariates: \n") 
              print(beta)
              cat("\nFrailty parameter: ", theta, "\n")
              if (length(object@knots) > 0) {
                  cat("\nInternal knots: \n") 
                  cat(object@knots, sep = ", ", fill = TRUE)
              }
              cat("\nBoundary knots: \n")
              cat(object@boundaryKnots, sep = ", ", fill = TRUE)
              if (object@degree > 0) {
                  cat("\nCoefficients of spline bases:\n")
                  print(alpha)
              } else {
                  cat("\nCoefficients of pieces:\n")
                  print(alpha)
              }
          })


#' @rdname show-method 
#' @aliases show,summaryRateReg-method
#' @importFrom stats printCoefmat
#' @export
setMethod(f = "show", signature = "summaryRateReg",
          definition = function(object) {
              if (attr(object@call, "show")) {
                  Call <- object@call
                  attr(Call, "show") <- NULL
                  cat("Call: \n")
                  print(Call)
                  cat("\n")
              }
              cat("Coefficients of covariates: \n") 
              printCoefmat(object@covarCoef)
              cat("\nParameter of frailty: \n")
              print(object@frailtyPar)
              if (attr(object@knots, "show")) {
                  if (length(object@knots) != 0) {
                      cat("\nInternal knots: \n")
                      cat(object@knots, sep = ", ", fill = TRUE)
                  }
                  cat("\nBoundary knots:\n")
                  cat(object@boundaryKnots, sep = ", ", fill = TRUE)
              }
              if (object@degree > 0) {
                  cat("\nDegree of spline bases:", object@degree, "\n")
                  cat("\nCoefficients of spline bases:\n")
                  printCoefmat(object@baseRateCoef)    
              } else {
                  cat("\nCoefficients of pieces:\n")
                  printCoefmat(object@baseRateCoef)    
              }
              cat("\nLoglikelihood: ", object@logL, "\n")
          })


#' @rdname show-method
#' @aliases show,sampleMcf-method 
#' @export
setMethod(f = "show", signature = "sampleMcf",
          definition = function(object) {
              cat("Call: \n")
              print(object@call)
              cat("\nFormula:\n")
              print(object@formula)
              cat("\nMCF:\n")
              print(object@MCF)
          })


#' @rdname show-method 
#' @aliases show,rateRegMcf-method
#' @export
setMethod(f = "show", signature = "rateRegMcf",
          definition = function(object) {
              cat("Formula:\n")
              print(object@formula)
              cat("\nNew data:\n")
              print(object@newdata)
              cat("\nConfidence level:",
                  paste(format(100 * object@level,
                               trim = TRUE, scientific = FALSE),
                        "%", sep = ""), "\n")
              cat("\nMCF:\n")
              print(object@MCF)
          })
