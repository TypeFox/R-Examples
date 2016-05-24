# @file Simulation.R
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

#' Simulate (negative) controls
#'
#' @details
#' Generate point estimates given known true effect sizes and standard errors
#'
#' @param n           Number of controls to simulate.
#' @param mean        The mean of the error distribution (on the log RR scale).
#' @param sd          The standard deviation of the error distribution (on the log RR scale).
#' @param seLogRr     The standard error of the log of the relative risk. This is recycled for the
#'                    controls. The default is to sample these from a uniform distribution.
#' @param trueLogRr   The true relative risk (on the log scale) used to generate these controls.  This
#'                    is recycled for the controls.
#'
#' @examples
#' data <- simulateControls(n = 50 * 3, mean = 0.25, sd = 0.25, trueLogRr = log(c(1, 2, 4)))
#' plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
#'
#' @export
simulateControls <- function(n = 50,
                             mean = 0,
                             sd = 0.1,
                             seLogRr = runif(n, min = 0.01, max = 0.2),
                             trueLogRr = 0) {
  theta <- rnorm(n, mean = mean, sd = sd)
  logRr <- rnorm(n, mean = trueLogRr + theta, sd = seLogRr)
  data <- data.frame(logRr = logRr, seLogRr = seLogRr, trueLogRr = trueLogRr)
  return(data)
}


# .simulation <- function() {
# 
#   negControls <- simulateControls()
#   null <- fitNull(negControls$logRr, negControls$seLogRr)
#   plotCalibrationEffect(negControls$logRr, negControls$seLogRr, null = null)
# 
#   data <- simulateControls(n = 50 * 3, trueLogRr = log(c(1, 2, 4)))
# 
#   plotTrueAndObserved(data$logRr, data$seLogRr, data$trueLogRr)
#   plotCoverage(data$logRr, data$seLogRr, data$trueLogRr)
# 
#   model <- fitSystematicErrorModel(data$logRr, data$seLogRr, data$trueLogRr)
# 
# 
#   cal <- calibrateConfidenceInterval(data$logRr, data$seLogRr, model)
# 
# 
#   plotTrueAndObserved(cal$logRr, cal$seLogRr, data$trueLogRr)
#   plotCoverage(cal$logRr, cal$seLogRr, data$trueLogRr)
# 
# 
#   logRr <- data$logRr
#   seLogRr <- data$seLogRr
#   trueLogRr <- data$trueLogRr
# 
#   data <- simulateControls()
#   data("sccs")
#   data <- sccs
#   p <- calibratePWithCiUsingMcmc(data$logRr,
#                                  data$seLogRr,
#                                  data$logRr[1],
#                                  data$seLogRr[1],
#                                  scale = c(0.05, 25),
#                                  iter = 10000)
#   mcmc <- attr(p, "mcmc")
#   mean(mcmc$acc)  # Acceptance rate
#   plot(ts(mcmc$chain[, 1]))  # Trace for the mean
#   plot(ts(mcmc$chain[, 2]))  # Trace for the precision (= 1/sqr(sd) )
#   mean(mcmc$chain[, 1])
#   mean(mcmc$chain[, 2])
#   1/(sqrt(quantile(mcmc$chain[, 2], c(0.025, 0.5, 0.975))))
#   p
#   null <- fitNull(data$logRr, data$seLogRr)
#   calibrateP(data$logRr[1], data$seLogRr[1], null, pValueConfidenceInterval = TRUE)
# }
