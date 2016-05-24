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


#' Akaike Information Criterion (AIC)
#' 
#' \code{AIC,rateReg-method} is an S4 class method calculating
#' Akaike information criterion (AIC) for one or several 
#' \code{rateReg-class} objects, according to the formula
#' - 2 * log-likelihood + 2 * nPar, where nPar represents the number
#' of parameters in the fitted model.
#'
#' When comparing models fitted by maximum likelihood to the same
#' data, the smaller the AIC, the better the fit.
#' \code{help(AIC, stats)} for other details.
#' 
#' @param object An object used to dispatch a method.
#' @param ... Optionally more fitted model objects.
#' @param k An optional numeric value used as the penalty per parameter.
#' The default \code{k = 2} is the classic AIC.
#' @return If just one object is provided, a numeric value representing
#' calculated AIC.
#' If multiple objects are provided, a data frame with rows
#' corresponding to the objects and columns \code{df} and \code{AIC},
#' where \code{df} means degree of freedom,
#' which is the number of parameters in the fitted model.
#' @aliases AIC,rateReg-method
#' @examples
#' ## See examples given in function rateReg.
#' @seealso
#' \code{\link{rateReg}} for model fitting;
#' \code{\link{summary,rateReg-method}} for summary of a fitted model;
#' \code{\link{BIC,rateReg-method}} for BIC.
#' @importFrom stats AIC
#' @export
setMethod(f = "AIC", signature = "rateReg",
          definition = function(object, ..., k = 2) {
              if (! missing(...)) {
                  inpList <- list(object, ...)

                  ## check on object class
                  checkRes <- sapply(inpList, objCheck)
                  if (any(! checkRes)) {
                      stop("Objects should be all from 'rateReg-class'.")
                  }

                  ## warning on different nObs
                  nObss <- sapply(inpList, nObsFun)
                  if (length(unique(nObss)) != 1) {
                      warning(paste("Models are not all fitted to the same",
                                    "number of observations."))
                  }
                  
                  abics <- sapply(inpList, abic, penal = k)
                  dfs <- sapply(inpList, sumDf)
                  val <- data.frame(df = dfs, AIC = abics)
                  Call <- match.call()
                  Call$k <- NULL
                  row.names(val) <- as.character(Call[-1L])
                  return(val)
              }
              ## else return
              abic(object, penal = k)
          })


#' Bayesian Information Criterion (BIC)
#' 
#' \code{BIC,rateReg-method} is an S4 class method calculating
#' Bayesian information criterion (BIC) or so-called
#' Schwarz's Bayesian criterion (SBC)
#' for one or several \code{rateReg-class} objects,
#' according to the formula
#' - 2 * log-likelihood + ln(nObs) * nPar,
#' where nPar represents the number of parameters in the fitted model
#' and nObs is the number of observations.
#' 
#' When comparing models fitted by maximum likelihood to the same
#' data, the smaller the BIC, the better the fit.
#' \code{help(BIC, stats)} for other details.
#' 
#' @param object An object used to dispatch a method.
#' @param ... Optionally more fitted model objects.
#' @return If just one object is provided, a numeric value representing
#' calculated BIC.
#' If multiple objects are provided, a data frame with rows
#' corresponding to the objects and columns \code{df} and \code{BIC},
#' where \code{df} means degree of freedom,
#' which is the number of parameters in the fitted model.
#' @aliases BIC,rateReg-method
#' @examples
#' ## See examples given in function rateReg.
#' @seealso
#' \code{\link{rateReg}} for model fitting;
#' \code{\link{summary,rateReg-method}} for summary of a fitted model;
#' \code{\link{AIC,rateReg-method}} for AIC.
#' @importFrom stats BIC
#' @export
setMethod(f = "BIC", signature = "rateReg",
          definition = function(object, ...) {
              if (! missing(...)) {
                  inpList <- list(object, ...)

                  ## check on object class
                  checkRes <- sapply(inpList, objCheck)
                  if (any(! checkRes)) {
                      stop("Objects should be all from 'rateReg-class'.")
                  }

                  nObss <- sapply(inpList, nObsFun)
                  k <- log(nObss)
                  abics <- sapply(seq_along(inpList), function (ind) {
                      abic(object = inpList[[ind]], penal = k[ind])
                  })
                  dfs <- sapply(inpList, sumDf)
                  val <- data.frame(df = dfs, BIC = abics)
                  Call <- match.call()
                  Call$k <- NULL
                  row.names(val) <- as.character(Call[-1L])
                  return(val)
              }
              ## else return
              abic(object, penal = log(object@nObs))
          })



### internal functions =========================================================
objCheck <- function (object) {
    inherits(object, "rateReg")
}

sumDf <- function(object) {
    sum(do.call("c", object@df))
}

abic <- function(object, penal) {
    - 2 * object@logL + penal * sumDf(object)
}

nObsFun <- function (object) {
    object@nObs
}

