# @file ConfidenceIntervalCalibration.R
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

#' Fit a systematic error model
#'
#' @details
#' Fit a model of the systematic error as a function of true effect size. This model is an extention
#' of the method for fitting the null distribution. The mean and standard deviations of the error
#' distributions are assumed to be linear with respect to the true effect size, and each component is
#' therefore represented by an intercept and a slope.
#'
#' @param logRr       A numeric vector of effect estimates on the log scale.
#' @param seLogRr     The standard error of the log of the effect estimates. Hint: often the standard
#'                    error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                    estimate>))/qnorm(0.025).
#' @param trueLogRr   A vector of the true effect sizes.
#'
#' @return
#' An object of type \code{systematicErrorModel}.
#'
#' @examples
#' controls <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(controls$logRr, controls$seLogRr, controls$trueLogRr)
#' model
#'
#' @export
fitSystematicErrorModel <- function(logRr, seLogRr, trueLogRr) {

  gaussianProduct <- function(mu1, mu2, sd1, sd2) {
    (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
  }

  LL <- function(theta, logRr, seLogRr, trueLogRr) {
    result <- 0
    for (i in 1:length(logRr)) {
      mean <- theta[1] + theta[2] * trueLogRr[i]
      sd <- theta[3] + theta[4] * trueLogRr[i]
      result <- result - log(gaussianProduct(logRr[i], mean, seLogRr[i], sd))
    }
    if (is.infinite(result))
      result <- 99999
    result
  }
  theta <- c(0, 1, 0.5, 0)
  fit <- optim(theta, LL, logRr = logRr, seLogRr = seLogRr, trueLogRr = trueLogRr, hessian = TRUE)
  fisher_info <- solve(fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  model <- fit$par
  names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  attr(model, "LB95CI") <- fit$par + qnorm(0.025) * prop_sigma
  attr(model, "UB95CI") <- fit$par + qnorm(0.975) * prop_sigma
  attr(model, "CovarianceMatrix") <- fisher_info
  class(model) <- "systematicErrorModel"
  model
}

#' Calibrate confidence intervals
#'
#' @details
#' Compute calibrated confidence intervals based on a model of the systematic error.
#'
#' @param logRr     A numeric vector of effect estimates on the log scale.
#' @param seLogRr   The standard error of the log of the effect estimates. Hint: often the standard
#'                  error = (log(<lower bound 95 percent confidence interval>) - log(<effect
#'                  estimate>))/qnorm(0.025).
#' @param model     An object of type \code{systematicErrorModel} as created by the
#'                  \code{\link{fitSystematicErrorModel}} function.
#'
#' @return
#' A data frame with calibrated confidence intervals and point estimates.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' model <- fitSystematicErrorModel(data$logRr, data$seLogRr, data$trueLogRr)
#' newData <- simulateControls(n = 15, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' result <- calibrateConfidenceInterval(newData$logRr, newData$seLogRr, model)
#' result
#'
#' @export
calibrateConfidenceInterval <- function(logRr, seLogRr, model) {
  statisticsConfidenceIntervals <- FALSE
  cal <- function(quantile, logRR, se, interceptLogRR, slopeLogRR, interceptSD, slopeSD) {
    z <- qnorm(quantile)
    (-sqrt((interceptLogRR - logRR)^2 * slopeSD^2 * z^2 - 2 * (interceptLogRR - logRR) * slopeLogRR *
      interceptSD * slopeSD * z^2 + slopeLogRR^2 * interceptSD^2 * z^2 + slopeLogRR^2 * se^2 *
      z^2 - slopeSD^2 * se^2 * z^4) - (interceptLogRR - logRR) * slopeLogRR + interceptSD * slopeSD *
      z^2)/(slopeLogRR^2 - slopeSD^2 * z^2)
  }
  result <- data.frame(logRr = rep(0, length(logRr)), logLb95Rr = 0, logUb95Rr = 0)
  for (i in 1:nrow(result)) {
    result$logRr[i] <- cal(0.5, logRr[i], seLogRr[i], model[1], model[2], model[3], model[4])
    result$logLb95Rr[i] <- cal(0.025, logRr[i], seLogRr[i], model[1], model[2], model[3], model[4])
    result$logUb95Rr[i] <- cal(0.975, logRr[i], seLogRr[i], model[1], model[2], model[3], model[4])
  }
  result$seLogRr <- (result$logLb95Rr - result$logRr)/qnorm(0.025)
  # Experimental: currently doesn't work because we get negative standard deviations
  if (statisticsConfidenceIntervals) {
    result$logRr_lb95 <- 0
    result$logRr_ub95 <- 0
    result$logLb95Rr_lb95 <- 0
    result$logLb95Rr_ub95 <- 0
    result$logUb95Rr_lb95 <- 0
    result$logUb95Rr_ub95 <- 0
    rand <- MASS::mvrnorm(10000, model, attr(model, "CovarianceMatrix"))
    for (i in 1:nrow(result)) {
      logRr <- cal(0.5, logRr[i], seLogRr[i], rand[1], rand[2], rand[3], rand[4])
      logLb95Rr <- cal(0.025, logRr[i], seLogRr[i], rand[1], rand[2], rand[3], rand[4])
      logUb95Rr <- cal(0.975, logRr[i], seLogRr[i], rand[1], rand[2], rand[3], rand[4])

      result$logRr_lb95 <- quantile(logRr, 0.025)
      result$logRr_ub95 <- quantile(logRr, 0.975)
      result$logLb95Rr_lb95 <- quantile(logLb95Rr, 0.025)
      result$logLb95Rr_ub95 <- quantile(logLb95Rr, 0.975)
      result$logUb95Rr_lb95 <- quantile(logUb95Rr, 0.025)
      result$logUb95Rr_ub95 <- quantile(logUb95Rr, 0.975)
    }
  }
  return(result)
}
