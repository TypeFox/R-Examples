# @file EmpiricalCalibrationUsingAsymptotics.R
#
# Copyright 2015 Observational Health Data Sciences and Informatics
#
# This file is part of EmpiricalCalibration
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Fit the null distribution
#'
#' @description
#' \code{fitNull} fits the null distribution to a set of negative controls
#'
#' @details
#' This function fits a Gaussian function to the negative control estimates as described in Schuemie
#' et al (2014).
#'
#'
#' @param logRr     A numeric vector of effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#'
#' @return
#' An object containing the parameters of the null distribution.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' null
#'
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#'
#' @export
fitNull <- function(logRr, seLogRr) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting null distribution")
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting null distribution")
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }

  gaussianProduct <- function(mu1, mu2, sd1, sd2) {
    (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
  }

  LL <- function(theta, estimate, se) {
    result <- 0
    for (i in 1:length(estimate)) {
      result <- result - log(gaussianProduct(estimate[i], theta[1], se[i], exp(theta[2])))
    }
    if (is.infinite(result))
      result <- 99999
    result
  }
  theta <- c(0, 0)
  fit <- optim(theta, LL, estimate = logRr, se = seLogRr)
  null <- fit$par
  null[2] <- exp(null[2])
  names(null) <- c("mean", "sd")
  class(null) <- "null"
  return(null)
}

#' @export
print.null <- function(x, ...) {
  writeLines("Estimated null distribution\n")
  output <- data.frame(Estimate = c(x[1], x[2]))
  colnames(output) <- c("Estimate")
  rownames(output) <- c("Mean", "SD")
  printCoefmat(output)

}

#' Calibrate the p-value
#'
#' @description
#' \code{calibrateP} computes calibrated p-values using the fitted null distribution
#'
#' @details
#' This function computes a calibrated two-sided p-value as described in Schuemie et al (2014).
#'
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#' @param null      An object of class \code{null} created using the \code{fitNull} function or an
#'                  object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.
#' @param ...       Any additional parameters (currently none).
#'
#' @return
#' The two-sided calibrated p-value.
#'
#' @examples
#' data(sccs)
#' negatives <- sccs[sccs$groundTruth == 0, ]
#' null <- fitNull(negatives$logRr, negatives$seLogRr)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' calibrateP(null, positive$logRr, positive$seLogRr)
#'
#' @references
#' Schuemie MJ, Ryan PB, Dumouchel W, Suchard MA, Madigan D. Interpreting observational studies: why
#' empirical calibration is needed to correct p-values. Statistics in Medicine 33(2):209-18,2014
#'
#' @export
calibrateP <- function(null, logRr, seLogRr, ...) {
  UseMethod("calibrateP")
}


#' @describeIn
#' calibrateP Computes the calibrated P-value using asymptotic assumptions.
#' @export
calibrateP.null <- function(null, logRr, seLogRr, ...) {

  oneAdjustedP <- function(logRR, se, null) {
    P_upper_bound <- pnorm((null[1] - logRR)/sqrt(null[2]^2 + se^2))
    P_lower_bound <- pnorm((logRR - null[1])/sqrt(null[2]^2 + se^2))
    2 * min(P_upper_bound, P_lower_bound)
  }

  adjustedP <- vector(length = length(logRr))
  for (i in 1:length(logRr)) {
    if (is.na(logRr[i]) || is.infinite(logRr[i]) || is.na(seLogRr[i]) || is.infinite(seLogRr[i])) {
      adjustedP[i] <- NA
    } else {
      adjustedP[i] <- oneAdjustedP(logRr[i], seLogRr[i], null)
    }
  }
  return(adjustedP)
}

#' Compute the (traditional) p-value
#'
#' @description
#' \code{computeTraditionalP} computes the traditional two-sided p-value based on the log of the
#' relative risk and the standerd error of the log of the relative risk.
#'
#' @param logRr     A numeric vector of one or more effect estimates on the log scale
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025)
#'
#' @return
#' The two-sided (traditional) p-value.
#'
#' @examples
#' data(sccs)
#' positive <- sccs[sccs$groundTruth == 1, ]
#' computeTraditionalP(positive$logRr, positive$seLogRr)
#'
#' @export
computeTraditionalP <- function(logRr, seLogRr) {
  z <- logRr/seLogRr
  return(2 * pmin(pnorm(z), 1 - pnorm(z)))  # 2-sided p-value
}
