################################################################################
################################################################################
## Crawford-Pendakur Upper and Lower Bound functions

# Copyright 2014 Julien Boelaert.
# 
# This file is part of revealedPrefs.
# 
# revealedPrefs is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# revealedPrefs is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with revealedPrefs.  If not, see <http://www.gnu.org/licenses/>.

## Function for Crawford-Pendakur-type Upper Bound algorithms
## for quantities x and prices p
cpUpper <- function(x, p, times= 1, afriat.par= 1, 
                    method= c("fastfloyd", "deep", "floyd")) {
  method <- match.arg(method)
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  the.samples <- replicate(times, sample(0:(nrow(x)-1)))
  
  if (method == "floyd") {
    the.call <- .Call("CpUp", p%*%t(x), the.samples, afriat.par, 
                      PACKAGE = "revealedPrefs")
  } else if (method == "fastfloyd") {
    the.call <- .Call("FastUp", p%*%t(x), the.samples, afriat.par,
                      PACKAGE = "revealedPrefs")
  } else 
    the.call <- .Call("DeepCpUp", x, p, the.samples, afriat.par,
                      PACKAGE = "revealedPrefs")
  
  the.call$clustering <- as.numeric(the.call$clustering)
  the.call$cluster.pop <- 
    as.numeric(the.call$cluster.pop[the.call$cluster.pop != 0])
  class(the.call) <- "upperBound"
  the.call$n.types <- length(the.call$cluster.pop)
  the.call$afriat.par <- afriat.par
  the.call
}

## Function for Crawford-Pendakur-type Lower Bound algorithms
## for quantities x and prices p
cpLower <- function(x, p, times= 1, afriat.par= 1) {
  if (any(is.na(x)) | any(is.na(p))) stop("NAs found in x or p\n")
  if (!all(dim(x) == dim(p))) stop("x and p must have same dimension\n")
  if (length(afriat.par) > 1 | afriat.par > 1 | afriat.par < 0)
    stop("'afriat.par' must be a real value between 0 and 1.\n")
  x <- as.matrix(x)
  p <- as.matrix(p)
  
  the.samples <- replicate(times, sample(0:(nrow(x)-1)))
  the.call <- .Call("CpLow", x, p, the.samples, afriat.par, 
                    PACKAGE= "revealedPrefs")
  the.call$violators <- the.call$violators + 1
  the.call$n.types <- length(the.call$violators)
  the.call$afriat.par <- afriat.par
  class(the.call) <- "lowerBound"
  the.call
}

################################################################################
################################################################################
## S3 methods

## Lower bound
print.lowerBound <- function(x, ...) {
  cat("  Lower bound on the number of types :", x$n.types, "\n")
}

summary.lowerBound <- function(object, ...) {
  cat("  Lower bound on the number of types :", object$n.types, "\n")
  cat("  (best of", length(object$hist.n.types), "run(s) of the algorithm)\n")
  cat("  Afriat parameter:", object$afriat.par,
      ifelse(object$afriat.par == 1, 
             "(no optimization error allowed)\n",
             paste("(", round(100 * (1 - object$afriat.par), 2), 
                   "% optimization error allowed)\n", sep= "")))
}

## Upper bound
print.upperBound <- function(x, ...) {
  cat("  Upper bound on the number of types :", x$n.types, "\n")
}

summary.upperBound <- function(object, ...) {
  cat("  Upper bound on the number of types :", object$n.types, "\n")
  cat("  Cluster populations                :", object$cluster.pop, "\n")
  cat("  (best of", length(object$hist.n.types), "run(s) of the algorithm)\n")
  cat("  Afriat parameter:", object$afriat.par,
      ifelse(object$afriat.par == 1, 
             "(no optimization error allowed)\n",
             paste("(", round(100 * (1 - object$afriat.par), 2), 
                   "% optimization error allowed)\n", sep= "")))
}
