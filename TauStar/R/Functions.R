# Copyright (C) 2015 Luca Weihs
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Computing t*
#'
#' Computes the t* U-statistic for input data pairs
#' (x_1,y_1), (x_2,y_2), ..., (x_n,y_n)
#' using the algorithm developed by Weihs, Drton, and Leung (2015)
#' <DOI:10.1007/s00180-015-0639-x>.
#'
#' @export
#'
#' @param x A numeric vector of x values (length >= 4).
#' @param y A numeric vector of y values, should be of the same length as x.
#' @param vStatistic If TRUE then will compute the V-statistic version of t*,
#'        otherwise will compute the U-Statistic version of t*. Default is to compute
#'        the U-statistic.
#' @param resample If TRUE then will compute an approximation of t* using a
#'        subsettting approach: samples of size sampleSize are taken from the data
#'        numResample times, t* is computed on each subsample, and all subsample t*
#'        values are then averaged. Note that this only works when vStatistic ==
#'        FALSE, in general you probably don't want to compute the V-statistic via
#'        resampling as the size of the bias depends on the sampleSize irrespective
#'        numResamples. Default is resample == FALSE so that t* is computed on all of
#'        the data, this may be slow for very large sample sizes.
#' @param numResamples See resample variable description for details, this
#'        value is ignored if resample == FALSE (ignored by default).
#' @param sampleSize See resample variable description for details, this value
#'        is ignored if resample == FALSE (ignored by default).
#' @param slow If TRUE then will compute the t* statistic using a naive O(n^4)
#'        algorithm. This option exists to allow for comparisons of the efficient new
#'        algorithm against the old slow one. At the moment this is not implemented to
#'        work when resample == TRUE.
#'
#' @return The numeric value of the t* statistic.
#'
#' @references
#' Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based
#' on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014),
#' no. 2, 1006--1028.
#' \cr\cr
#' Weihs, Luca, Mathias Drton, and Dennis Leung. "Efficient Computation of the
#' Bergsma-Dassios Sign Covariance." arXiv preprint arXiv:1504.00964 (2015).
#'
#' @examples
#' \dontrun{
#' library(TauStar)
#'
#' # Compute t* for a concordant quadruple
#' tStar(c(1,2,3,4), c(1,2,3,4)) # == 2/3
#'
#' # Compute t* for a discordant quadruple
#' tStar(c(1,2,3,4), c(1,-1,1,-1)) # == -1/3
#'
#' # Compute t* on random normal iid normal data
#' set.seed(23421)
#' tStar(rnorm(4000), rnorm(4000)) # near 0
#'
#' # Compute t* as a v-statistic
#' set.seed(923)
#' tStar(rnorm(100), rnorm(100), vStatistic=TRUE)
#'
#' # Compute an approximation of tau* via resampling
#' set.seed(9492)
#' tStar(rnorm(10000), rnorm(10000), resample=TRUE, sampleSize=30,
#'       numResamples=5000)
#' }
tStar <- function(x, y, vStatistic = FALSE, resample = FALSE,
                  numResamples = 500, sampleSize = min(length(x), 1000),
                  slow = FALSE) {
  if(!is.numeric(x) || !is.numeric(y)) {
    stop("Input x and y to tStar must be numeric.")
  }
  if(length(x) != length(y) || length(x) < 4) {
    stop("Input x and y to tStar are of the wrong length, they must both have equal length < 4.")
  }
  if(!is.logical(vStatistic) || length(vStatistic) != 1) {
    stop("Input parameter vStatistic into function tStar must be a logical T/F value.")
  }
  if(!is.logical(slow) || length(slow) != 1) {
    stop("Input parameter slow into function tStar must be a logical T/F value.")
  }
  if(!is.logical(resample) || length(resample) != 1) {
    stop("Input parameter resample into function tStar must be a logical T/F value.")
  }
  if(resample && numResamples <= 0) {
    stop("When resampling the number of resamples must be positive.")
  }
  if(resample && (sampleSize < 4 || sampleSize > length(x))) {
    stop("When resampling the sample size must be greater than 3 and less than the length of the input data.")
  }
  if(resample && slow) {
    stop("Resampling is not currently implemented with the slow algorithm.")
  }
  if(resample && vStatistic) {
    stop("Resampling is not currently implemented when computing V-statistics. Note that you probably don't want to compute the V-statistic via resampling as the size of the bias would depend on the size of subsets chosen independent of the number of resamples.")
  }
  ord = sort.list(x, method = "quick", na.last = NA)
  x = x[ord]
  y = y[ord]
  if(resample) {
    return(TStarFastResampleRCPP(x, y, numResamples, sampleSize))
  } else if(slow) {
    return(TStarSlowTiesRCPP(x, y, vStatistic))
  } else if (vStatistic) {
    return(VTStarFastTiesRCPP(x, y))
  } else {
    return(TStarFastTiesRCPP(x, y))
  }
}

#' Quantiles of a distribution.
#'
#' Computes the pth quantile of a cumulative distribution function using a
#' simple binary serach algorithm. This can be extremely slow but has the
#' benefit of being trivial to implement.
#'
#' @param pDistFunc a cumulative distribution function on the real numbers, it
#'        should take a single argument x and return the cumualtive distribution
#'        function evaluated at x.
#' @param p the quantile \eqn{p\in[0,1]}
#' @param lastLeft binary search works by continuously decreasing the search
#'        space from the left and right. lastLeft should be a lower bound for
#'        the quantile p.
#' @param lastRight similar to lastRight but should be an upper bound.
#' @param error the error tolerated from the binary search
#'
#' @return the quantile (within error).
binaryQuantileSearch = function(pDistFunc, p, lastLeft, lastRight,
                                error = 10^-4) {
  mid = (lastLeft + lastRight) / 2
  pMid = pDistFunc(mid)
  if (lastRight - lastLeft < error) {
    return(mid)
  } else if (pMid > p) {
    return(binaryQuantileSearch(pDistFunc, p, lastLeft, mid, error))
  } else {
    return(binaryQuantileSearch(pDistFunc, p, mid, lastRight, error))
  }
}

#' Null asymptotic distribution of t* in the continuous case
#'
#' Density, distribution function, quantile function and random generation for
#' the asymptotic null distribution of t* in the continuous case. That is, in
#' the case that t* is generated from a sample of jointly continuous independent
#' random variables.
#'
#' @export
#'
#' @param x the value (or vector of values) at which to evaluate the function.
#' @param p the probability (or vector of probabilities) for which to get the
#'        quantile.
#' @param n the number of observations to return.
#' @param error a tolerated error in the result. This should be considered as a
#'        guide rather than an exact upper bound to the amount of error.
#' @param lower.tail a logical value, if TRUE (default), probabilities are
#'        \eqn{P(X\leq x)} otherwise \eqn{P(X>x)}.
#'
#' @name pHoeffInd
#' @rdname pHoeffInd
#'
#' @return dHoeffInd gives the density, pHoeffInd gives the distribution
#'         function, qHoeffInd gives the quantile function, and rHoeffInd
#'         generates random samples.
pHoeffInd <- function(x, lower.tail = T, error = 10^-5) {
  if (lower.tail) {
    return(HoeffIndCdfRCPP(x + 1, error))
  }
  return(1 - HoeffIndCdfRCPP(x + 1, error))
}

#' @rdname pHoeffInd
#' @export
rHoeffInd <- function(n) {
  sims = numeric(n)
  for(i in 1:50) {
    for(j in 1:50) {
      sims = sims + (36 / pi^4) * (1 / (i^2 * j^2)) * (rchisq(n, df = 1) - 1)
    }
  }
  return(sims)
}

#' @rdname pHoeffInd
#' @export
dHoeffInd <- function(x, error = 1/2 * 10^-3) {
  if (length(x) != 1) {
    return(sapply(x, function(y) { HoeffIndPdfRCPP(y + 1, error) }))
  }
  HoeffIndPdfRCPP(x + 1, error)
}

#' @rdname pHoeffInd
#' @export
qHoeffInd <- function(p, error = 10^-4) {
  if (p == 0) {
    return(-1)
  } else if (p == 1) {
    return(Inf)
  } else if (p < 0 || p > 1) {
    stop("p must be between 0 and 1 in qHoeffInd.")
  }
  right = 1
  while (pHoeffInd(right, error = error / 2) < p) {
    right = right * 2
  }
  pDistFunc = function(x) { pHoeffInd(x, error/2) }
  return(binaryQuantileSearch(pDistFunc, p, -1,
                              right, error / 2))
}

#' Check if Vector of Probabilities
#'
#' Checks if the input vector has entries that sum to 1 and are non-negative
#'
#' @param probs the probability vector to check
#'
#' @return TRUE if conditions are met, FALSE if otherwise
isProbVector <- function(probs) {
  probSum = sum(probs)
  return(abs(probSum - 1) < 10^-7 && all(probs >= 0))
}

#' Check if a Valid Probability
#'
#' Checks if the input vector has a single entry that is between 0 and 1
#'
#' @param prob the probability to check
#'
#' @return TRUE if conditions are met, FALSE if otherwise
isProb <- function(prob) {
  return(length(prob) == 1 && prob >= 0 && prob <= 1)
}

#' Null asymptotic distribution of t* in the discrete case
#'
#' Density, distribution function, quantile function and random generation for
#' the asymptotic null distribution of t* in the discrete case. That is, in the
#' case that t* is generated from a sample of jointly discrete independent
#' random variables X and Y.
#'
#' @export
#'
#' @inheritParams pHoeffInd
#'
#' @param probs1 a vector of probabilities corresponding to the (ordered)
#'        support of X. That is if your first random variable has support
#'        \eqn{u_1,...,u_n} then the ith entry of probs should be
#'        eqn{P(X = u_i)}.
#' @param probs2 just as probs1 but for the second random variable Y.
#'
#' @name pDisHoeffInd
#' @rdname pDisHoeffInd
#'
#' @return dDisHoeffInd gives the density, pDisHoeffInd gives the distribution
#'         function, qDisHoeffInd gives the quantile function, and
#'         rDisHoeffInd generates random samples.
pDisHoeffInd <- function(x, probs1, probs2, lower.tail = T, error = 10^-5) {
  if (!isProbVector(probs1) || !isProbVector(probs2)) {
    stop("either probs1 or probs2 in pDisHoeffInd is not a probability vector.")
  }
  if (any(probs1 == 1) || any(probs2 == 1)) {
    lowerTailProb = if (x >= 0) 1 else 0
  } else {
    eigenP = eigenForDiscreteProbs(probs1)
    eigenQ = eigenForDiscreteProbs(probs2)
    lowerTailProb = HoeffIndDiscreteCdfRCPP(x + 4 * sum(eigenP) * sum(eigenQ),
                                            eigenP, eigenQ, error)
  }
  if (lower.tail) {
    return(lowerTailProb)
  }
  return(1 - lowerTailProb)
}

#' @rdname pDisHoeffInd
#' @export
dDisHoeffInd <- function(x, probs1, probs2, error = 10^-3) {
  if (!isProbVector(probs1) || !isProbVector(probs2)) {
    stop("either probs1 or probs2 in pDisHoeffInd is not a probability vector.")
  }
  if (any(probs1 == 1) || any(probs2 ==1)) {
    stop(paste("since length(probs1) == 1 or length(probs2) == 1 resulting",
         "distribution is degenerate and thus has no density."))
  }
  eigenP = eigenForDiscreteProbs(probs1)
  eigenQ = eigenForDiscreteProbs(probs2)
  return(HoeffIndDiscretePdfRCPP(x + 4 * sum(eigenP) * sum(eigenQ),
                                 eigenP,
                                 eigenQ,
                                 error))
}

#' @rdname pDisHoeffInd
#' @export
rDisHoeffInd <- function(n, probs1, probs2) {
  if (!isProbVector(probs1) || !isProbVector(probs2)) {
    stop("either probs1 or probs2 in rDisHoeffInd is not a probability vector.")
  }
  if (any(probs1 == 1) || any(probs2 ==1)) {
    return(rep(0, n))
  }
  eigenP = eigenForDiscreteProbs(probs1)
  eigenQ = eigenForDiscreteProbs(probs2)
  asymResults = numeric(n)
  for (i in 1:length(eigenP)) {
    for (j in 1:length(eigenQ)) {
      asymResults = asymResults + 4 * eigenP[i] * eigenQ[j] * (rchisq(n, df = 1) - 1)
    }
  }
  return(asymResults)
}

#' @rdname pDisHoeffInd
#' @export
qDisHoeffInd <- function(p, probs1, probs2, error = 10^-4) {
  if (!isProbVector(probs1) || !isProbVector(probs2)) {
    stop("either probs1 or probs2 in pqisHoeffInd is not a probability vector.")
  }
  if (!isProb(p)) {
    stop("probability p must be length 1 and in [0,1].")
  }
  if (any(probs1 == 1) || any(probs2 ==1)) {
    return(0)
  }

  left = -4 * sum(eigenForDiscreteProbs(probs1)) * sum(eigenForDiscreteProbs(probs2))
  if (p == 0) {
    return(left)
  } else if (p == 1) {
    return(Inf)
  }
  right = 1
  while (pDisHoeffInd(right, probs1, probs2, error = error / 2) < p) {
    right = right * 2
  }
  return(binaryQuantileSearch(function(x) { pDisHoeffInd(x, probs1, probs2,
                                                         error = error / 2) },
                              p, left, right, error))
}

#' Null asymptotic distribution of t* in the mixed case
#'
#' Density, distribution function, quantile function and random generation for
#' the asymptotic null distribution of t* in the mixed case. That is, in the
#' case that t* is generated a sample from an independent bivariate distribution
#' where one coordinate is marginally discrete and the other marginally
#' continuous.
#'
#' @export
#'
#' @inheritParams pHoeffInd
#'
#' @param probs a vector of probabilities corresponding to the (ordered)
#'        support the marginally discrete random variable. That is, if the
#'        marginally discrete distribution has support \eqn{u_1,...,u_n}
#'        then the ith entry of probs should be the probability of seeing
#'        \eqn{u_i}.
#'
#' @name pMixHoeffInd
#' @rdname pMixHoeffInd
#'
#' @return dMixHoeffInd gives the density, pMixHoeffInd gives the distribution
#'         function, qMixHoeffInd gives the quantile function, and
#'         rMixHoeffInd generates random samples.
pMixHoeffInd <- function(x, probs, lower.tail = T, error = 10^-6) {
  if (!isProbVector(probs)) {
    stop("probs in pMixHoeffInd is not a probability vector.")
  }
  if (any(probs == 1)) {
    lowerTailProb = if (x >= 0) 1 else 0
  } else {
    eigenP = eigenForDiscreteProbs(probs)
    lowerTailProb = HoeffIndMixedCdfRCPP(x + 2 * sum(eigenP), eigenP, error)
  }
  if (lower.tail) {
    return(lowerTailProb)
  }
  return(1 - lowerTailProb)
}

#' @rdname pMixHoeffInd
#' @export
dMixHoeffInd <- function(x, probs, error = 10^-3) {
  if (!isProbVector(probs)) {
    stop("probs in dMixHoeffInd is not a probability vector.")
  }
  if (any(probs == 1)) {
    stop(paste("probs represents a degenerate discrete random variable and",
               "thus the density of the asymptotic distribution doesn't exist"))
  }
  eigenP = eigenForDiscreteProbs(probs)
  return(HoeffIndMixedPdfRCPP(x + 2 * sum(eigenP),
                              eigenP,
                              error))
}

#' @rdname pMixHoeffInd
#' @export
rMixHoeffInd <- function(n, probs, error = 10^-8) {
  if (!isProbVector(probs)) {
    stop("probs in rMixHoeffInd is not a probability vector.")
  }
  if (any(probs == 1)) {
    return(rep(0, n))
  }
  eigenP = eigenForDiscreteProbs(probs)
  top = ceiling((sum((12 / pi^2 * eigenP)^2) / (3 * error))^(1/3) + 1)
  sims = numeric(n)
  for (lambda in eigenP) {
    for (i in 1:top) {
      sims = sims + (12 / pi^2) * lambda/i^2 * (rchisq(n, df = 1) - 1)
    }
  }
  return(sims)
}

#' @rdname pMixHoeffInd
#' @export
qMixHoeffInd <- function(p, probs, error = 10^-4) {
  if (!isProbVector(probs)) {
    stop("probs in dMixHoeffInd is not a probability vector.")
  }
  if (!isProb(p)) {
    stop("probability p must be length 1 and in [0,1].")
  }
  if (any(probs == 1)) {
    return(0)
  }
  left = -2 * sum(eigenForDiscreteProbs(probs))
  if (p == 0) {
    return(left)
  } else if (p == 1) {
    return(Inf)
  }
  right = 1
  while (pMixHoeffInd(right, probs, error = error / 2) < p) {
    right = right * 2
  }
  return(binaryQuantileSearch(function(x) { pMixHoeffInd(x, probs,
                                                         error = error / 2) },
                              p, left, right, error))
}

#' Print Tau* Test Results
#'
#' A simple print function for tstest (Tau* test) objects.
#'
#' @param tsObj the tstest object to be printed
print.tstest <- function(tsObj) {
  asymTest = F
  if (tsObj$mode %in% c("continuous", "discrete", "mixed")) {
    cat(paste("Test Type: asymptotic", tsObj$mode, "\n"))
    asymTest = T
  } else {
    cat(paste("Test Type: permutation (", tsObj$resamples," simulations)\n",
              sep = ""))
  }
  cat(paste("Input Length:", length(tsObj$x), "\n\n"))

  cat(paste("Results:\n"))
  df = data.frame(round(tsObj$tStar, 5))
  if (asymTest) {
    df = cbind(df, round(tsObj$pVal, 5))
    colnames(df) = c("t* value", "Asym. p-val")
  } else {
    df = cbind(df, round(tsObj$permPVal, 5))
    colnames(df) = c("t* value", "Perm. p-val")
  }
  row.names(df) = ""
  print(df)
}

#' Determine if input data is discrete
#'
#' Attempts to determine if the input data is from a discrete distribution. Will
#' return true if the data type is of type integer or there are non-unique
#' values.
#'
#' @param x a vector which should be determined if discrete or not.
#'
#' @return the best judgement of whether or not the data was discrete
isDiscrete = function(x) {
  if (is.integer(x) || (length(unique(x)) != length(x))) {
    return(T)
  }
  return(F)
}

#' Is Vector Valid Data?
#'
#' Determines if input vector is a valid vector of real valued observations
#'
#' @param x the vector to be tested
#'
#' @return TRUE or FALSE
isValidDataVector <- function(x) {
  return(is.numeric(x) || is.integer(x))
}

#' Test of Independence Using the Tau* Measure
#'
#' Performs a (consistent) test of independence between two input vectors using
#' the asymptotic (or permutation based) distribution of the test statistic t*.
#' The asymptotic results hold in the case that x is generated from either a
#' discrete or continous distribution and similarly for y (in particular it is
#' allowed for one to be continuous while the other is discrete). The asymptotic
#' distributions were computed in Nandy, Weihs, and Drton (2016)
#' <http://arxiv.org/abs/1602.04387>.
#'
#' @export
#'
#' @param x a vector of sampled values.
#' @param y a vector of sampled values corresponding to x, y must be the same
#'        length as x.
#' @param mode should be one of five possible values: "auto", "continuous",
#'        "discrete", "mixed", or "permutation". If "auto" is selected then the
#'        function will attempt to automatically determine whether x,y are
#'        discrete or continuous and then perform the appropriate asymptotic
#'        test. In cases "continuous", "discrete", and "mixed" we perform the
#'        associated asymptotic test making the given assumption. Finally
#'        if "permutation" is selected then the function runs a Monte-Carlo
#'        permutation test for some given number of resamplings.
#' @param resamples the number of resamplings to do if mode = "permutation".
#'        Otherwise this value is ignored.
#'
#' @return a list with class "tstest" recording the outcome of the test.
#'
#' @references
#' Preetam Nandy, Luca Weihs, and Mathias Drton. Large-Sample Theory for the
#' Bergsma-Dassios Sign Covariance. arXiv preprint arXiv:1602.04387. 2016.
#'
#' @examples
#' set.seed(123)
#' x = rnorm(100)
#' y = rnorm(100)
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # big p-value
#'
#' y = y + x # make x and y correlated
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # small p-value
tauStarTest <- function(x, y, mode="auto", resamples = 1000) {
  if (!isValidDataVector(x) || !isValidDataVector(y) ||
      length(x) != length(y)) {
    stop(paste("vectors inputted to tauStarTest must be of type numeric or",
               "integer and must be the same length"))
  }
  if (length(resamples) != 1 || resamples %% 1 != 0) {
    stop("resamples must be integer valued.")
  }
  xIsDis = isDiscrete(x)
  yIsDis = isDiscrete(y)
  n = length(x)

  toReturn = list()
  class(toReturn) = "tstest"
  toReturn$mode = mode
  toReturn$x = x
  toReturn$y = y
  toReturn$tStar = tStar(x, y)
  toReturn$resamples = resamples

  if (mode == "auto") {
    if (xIsDis && yIsDis) {
      mode = "discrete"
    } else if (xIsDis || yIsDis) {
      mode = "mixed"
    } else {
      mode = "continuous"
    }
    toReturn$mode = paste(toReturn$mode, mode, sep = "-")
  }

  if (mode == "continuous") {
    if (xIsDis || yIsDis) {
      stop("Input vectors to tauStarTest are have repeated entries or are of the integer type but the mode is set to continuous.")
    }
    toReturn$pVal = 1 - pHoeffInd(n * toReturn$tStar)

  } else if (mode == "discrete") {
    p = as.numeric(table(x)) / n
    q = as.numeric(table(y)) / n
    toReturn$pVal = 1 - pDisHoeffInd(n * toReturn$tStar, probs1 = p, probs2 = q)

  } else if (mode == "mixed") {
    if (xIsDis && yIsDis) {
      stop("Input vectors to tauStarTest both are discrete but 'mixed' mode was chosen instead of 'discrete'.\n")
    } else if (!xIsDis && !yIsDis) {
      warning(paste("Neither vector input to tauStarTest has duplicate entries",
                    "but 'mixed' mode was selected. Will default to assuming x",
                    "is discrete and y continuous. \n"))
    }
    if (xIsDis) {
      z = x
      x = y
      y = z
    }
    p = as.numeric(table(y)) / n
    toReturn$pVal = 1 - pMixHoeffInd(n * toReturn$tStar, probs = p)
  } else if (mode == "permutation") {
    sampleTStars = numeric(resamples)
    for (i in 1:resamples) {
      sampleTStars[i] = tStar(sample(x), y)
    }
    toReturn$pVal = mean(sampleTStars >= toReturn$tStar)
  } else {
    stop("Invalid mode as input to tauStarTest.")
  }
  return(toReturn)
}
