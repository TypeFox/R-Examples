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


#' Summarizing a Fitted Model
#'
#' Summary of estimated coefficients of covariates, rate function bases,
#' and estimated rate parameter of frailty random variable, etc.,
#' which can be printed out by \code{show}.
#' 
#' \code{summary,rateReg-method} returns a
#' \code{\link{summaryRateReg-class}} object,
#' whose slots include
#' \itemize{
#'     \item \code{covarCoef}: Estimated covariate coefficients.
#'     \item \code{frailtyPar}: Estimated rate parameter of gamma frailty.
#'     \item \code{baseRateCoef}: Estimated coeffcients of baseline
#'         rate function.
#' }
#' For the meaning of other slots, see \code{\link{rateReg}}.
#' 
#' @param object \code{\link{rateReg-class}} object.
#' @param showCall A logic value with dafault \code{TRUE},
#' indicating whether function \code{show} 
#' prints out the original call information of \code{rateReg}.
#' It may be helpful for a more concise printout.
#' @param showKnots A logic value with default \code{TRUE}, 
#' indicating whether function \code{show}
#' prints out the internal and boundary knots.
#' Similar to argument \code{showCall}, It may be helpful
#' for a more concise printout.
#' @param ... Other arguments for future usage.
#' @return summaryRateReg-class object
#' @aliases summary,rateReg-method
#' @examples
#' ## See examples given in function rateReg.
#' @seealso \code{\link{rateReg}} for model fitting;
#' \code{\link{coef,rateReg-method}} for point estimates of
#' covariate coefficients; 
#' \code{\link{confint,rateReg-method}} for confidence intervals
#' of covariate coeffcients;
#' \code{\link{baseRate,rateReg-method}} for coefficients of baseline
#' rate function.
#' @export
setMethod(f = "summary", signature = "rateReg",
          definition = function(object, showCall = TRUE,
                                showKnots = TRUE, ...) {
              Call <- object@call
              attr(Call, "show") <- showCall
              knots <- object@knots
              boundaryKnots <- object@boundaryKnots
              attr(knots, "show") <- showKnots
              beta <- object@estimates$beta
              theta <- object@estimates$theta
              alpha <- object@estimates$alpha
              ## check on object validity by 'new', validObject(results)
              results <- new("summaryRateReg", 
                             call = Call,
                             knots = knots,
                             boundaryKnots = boundaryKnots,
                             covarCoef = beta,
                             frailtyPar = theta,
                             degree = object@degree,
                             baseRateCoef = alpha,
                             logL = object@logL)
              ## return
              results
          })

