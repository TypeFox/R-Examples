# @file EmpiricalCalibrationUsingPLikelihood.R
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

# This doesn't work in all cases, so abandoning it for now
calibratePUsingPLikelihood <- function(logRrNegatives,
                                       seLogRrNegatives,
                                       logRrPositives,
                                       seLogRrPositives) {
  gaussianProduct <- function(mu1, mu2, sd1, sd2) {
    (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
  }

  logLikelihoodMuSigma <- function(mu, sigma, logRr, seLogRr) {
    result <- 0
    for (i in 1:length(logRr)) {
      result <- result + log(gaussianProduct(logRr[i], mu, seLogRr[i], sigma))
    }
    if (is.infinite(result)) {
      return(-99999)
    } else {
      return(result)
    }
  }

  binarySearchSigma <- function(p, mu, logRr, seLogRr, upperBound = TRUE, precision = 1e-07) {
    logP <- log(p)
    if (logRr > mu) {
      flip <- TRUE
    } else {
      flip <- FALSE
    }
    if (!upperBound) {
      flip <- !flip
    }
    L <- -1000
    H <- 2

    while (H >= L) {
      M <- L + (H - L)/2
      # print(M)
      if (upperBound) {
        target <- pnorm((mu - logRr)/sqrt(exp(M)^2 + seLogRr^2), log.p = TRUE) - logP
      } else {
        target <- pnorm((logRr - mu)/sqrt(exp(M)^2 + seLogRr^2), log.p = TRUE) - logP
      }
      # print(target)
      if (target > precision) {
        if (flip) {
          H <- M
        } else {
          L <- M
        }
      } else if (-target > precision) {
        if (flip) {
          L <- M
        } else {
          H <- M
        }
      } else {
        return(exp(M))
      }
      if (M == -1000 || M == 2)
        return(0)
    }
  }

  likelihoodMu <- function(mu, p, logRr, seLogRr, logRrNegatives, seLogRrNegatives, upperBound) {
    lik <- vector(length = length(mu))
    for (i in 1:length(mu)) {
      # writeLines(paste('mu',mu[i]))
      sigma <- binarySearchSigma(p, mu[i], logRr, seLogRr, upperBound)
      # writeLines(paste('sigma',sigma))
      if (sigma == 0) {
        lik[i] <- 0
      } else {
        lik[i] <- exp(logLikelihoodMuSigma(mu[i], sigma, logRrNegatives, seLogRrNegatives))
      }
      # writeLines(paste('lik',lik[i]))
    }
    return(lik)
  }

  logLikelihoodP <- function(p,
                             logRr,
                             seLogRr,
                             logRrNegatives,
                             seLogRrNegatives,
                             upperBound = TRUE) {
    lik <- integrate(likelihoodMu,
                     lower = -10,
                     upper = 10,
                     p = p,
                     logRr = logRr,
                     seLogRr = seLogRr,
                     logRrNegatives = logRrNegatives,
                     seLogRrNegatives = seLogRrNegatives,
                     upperBound = upperBound,
                     subdivisions = 1000)$value
    if (lik == 0) {
      print(paste(p, -9999999))
      return(-9999999)
    } else {
      print(paste(p, log(lik)))
      return(log(lik))
    }
  }

  profile <- function(mode,
                      alpha = 0.025,
                      upperBound = FALSE,
                      logRr = logRr,
                      seLogRr = seLogRr,
                      logRrNegatives,
                      seLogRrNegatives,
                      precision = 1e-04,
                      minP = 1e-06) {
    q <- qchisq(1 - alpha, 1)/2
    if (!upperBound) {
      # flip <- FALSE
      L <- 0
      H <- mode
      q <- qchisq(1 - alpha, 1)/2
    } else {
      # flip <- TRUE
      L <- mode
      H <- 0.5
      q <- qchisq(alpha, 1)/2
    }
    llMode <- logLikelihoodP(mode,
                             logRr = logRr,
                             seLogRr = seLogRr,
                             logRrNegatives = logRrNegatives,
                             seLogRrNegatives = seLogRrNegatives)

    # Determine whether gradient goes up or down:
    llModeDelta <- logLikelihoodP(mode * 0.1,
                                  logRr = logRr,
                                  seLogRr = seLogRr,
                                  logRrNegatives = logRrNegatives,
                                  seLogRrNegatives = seLogRrNegatives)
    if (llModeDelta > llMode) {
      flip <- !upperBound
    } else {
      flip <- upperBound
    }

    previousMetric <- 0
    while (H >= L) {
      M <- L + (H - L)/2
      llM <- logLikelihoodP(M,
                            logRr = logRr,
                            seLogRr = seLogRr,
                            logRrNegatives = logRrNegatives,
                            seLogRrNegatives = seLogRrNegatives)
      metric <- llMode - llM - q
      if (metric == previousMetric) {
        # Reached precision of the likelihood function. Assume we're already close:
        return(L)
      }
      previousMetric <- metric
      # writeLines(paste('M =', M, 'Metric = ',metric))
      if (metric > precision) {
        if (flip) {
          H <- M
        } else {
          L <- M
        }
      } else if (-metric > precision) {
        if (flip) {
          L <- M
        } else {
          H <- M
        }
      } else {
        return(M)
      }
      if (M < minP) {
        # Could be a long way down, but we don't really care
        return(0)
      }
    }
    return(L)
  }

  adjustedP <- data.frame(p = rep(1, length(logRrPositives)), lb95ci = 0, ub95ci = 0)
  for (i in 1:length(logRrPositives)) {
    pUpperBound <- optimize(logLikelihoodP,
                            c(0, 1),
                            upperBound = TRUE,
                            logRr = logRrPositives[i],
                            seLogRr = seLogRrPositives[i],
                            logRrNegatives = logRrNegatives,
                            seLogRrNegatives = seLogRrNegatives,
                            maximum = TRUE,
                            tol = 1e-04)$maximum
    pLowerBound <- optimize(logLikelihoodP,
                            c(0, 1),
                            upperBound = FALSE,
                            logRr = logRrPositives[i],
                            seLogRr = seLogRrPositives[i],
                            logRrNegatives = logRrNegatives,
                            seLogRrNegatives = seLogRrNegatives,
                            maximum = TRUE,
                            tol = 1e-04)$maximum
    p <- min(pUpperBound, pLowerBound)
    upperBound <- profile(p,
                          alpha = 0.025,
                          upper = TRUE,
                          logRr = logRrPositives[i],
                          seLogRr = seLogRrPositives[i],
                          logRrNegatives = logRrNegatives,
                          seLogRrNegatives = seLogRrNegatives)
    lowerBound <- profile(p,
                          alpha = 0.025,
                          upper = FALSE,
                          logRr = logRrPositives[i],
                          seLogRr = seLogRrPositives[i],
                          logRrNegatives = logRrNegatives,
                          seLogRrNegatives = seLogRrNegatives)
    adjustedP$p[i] <- p * 2
    adjustedP$lb95ci[i] <- lowerBound * 2
    adjustedP$ub95ci[i] <- upperBound * 2
  }
  return(adjustedP)
}
