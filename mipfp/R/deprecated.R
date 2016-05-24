# File mipfp/R/deprecated.R
# by Johan Barthelemy and Thomas Suesse
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ------------------------------------------------------------------------------
# This file contains the deprecated methods that will be removed in future
# releases of the package.
# ------------------------------------------------------------------------------

IpfpCov <- function(estimate, seed, target.list, replace.zeros = 1e-10) {
  # Compute the covariance matrix of the estimators produced by Ipfp.
  #
  # This function determines the covariance matrix of the estimated proportions
  # using the Delta method given in the paper "Models for Contingency Tables 
  # With Known Margins When Target and Sampled Populations Differ" written by 
  # Little and Wu (1991).
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   estimate: The array of estimate produced by the Ipfp function.
  #   seed: The initial multi-dimensional array updated by Ipfp. 
  #   target.list: A list of the target margins used by the Ipfp function. Each
  #                component of the list is an array whose cells indicates
  #                which dimension the corresponding margin relates to.
  #   replace.zeros: If 0-cells are to be found in either the seed or the
  #                  estimate arrays, then their values are replaced with this
  #                  value.
  #
  # Returns: A list whose elements are  
  #   p.hat.cov: A covariance matrix of the estimated probabilities (last index 
  #              move fastest).
  #   df: Degrees of freedom of the estimates.
  
  n <- sum(seed)  
  seed.prob <- Array2Vector(seed / sum(seed))
  estimate.prob <- Array2Vector(estimate / sum(estimate))
  
  # checking if 0-cells values and replace them with a small value  
  seed.prob <- ifelse(seed.prob == 0, replace.zeros, seed.prob)    
  estimate.prob <- ifelse(estimate.prob == 0, replace.zeros, estimate.prob)
  
  # computation of the diagonal matrix filled with the inverses of seed and
  # estimated probabilities
  D.estimate <- diag(1 / estimate.prob)
  D.seed <- diag(1 / seed.prob)
  
  # computation of A such that A' * vector(estimate) = vector(target.data)
  # ... one line filled with ones
  A.transp <- matrix(1, nrow = 1, ncol = length(estimate.prob))  
  # ... constrainst (removing the first one since it is redundant information)
  for (j in 1:length(target.list)) {
    marg.mat <- cmm::MarginalMatrix(var = 1:length(dim(seed)), 
                                    marg = target.list[[j]], 
                                    dim = dim(seed))[-1, ]
    
    A.transp <- rbind(marg.mat, A.transp)
  }  
  A <- t(A.transp)
  
  # removing the linearly dependant columns from A (redundant constrainst)
  A <- GetLinInd(A)$mat.li
  
  # computation of the orthogonal complement of A (using QR decomposition)
  K <- qr.Q(qr(A), complete = TRUE)[, (dim(A)[2] + 1):dim(A)[1]]
  
  # computation of the variance  
  estimate.var <- (1 / n) * K %*% solve((t(K) %*% D.estimate %*% K)) %*%
    t(K) %*% D.seed %*% K %*%
    solve(t(K) %*% D.estimate %*% K) %*% t(K)   
  
  # returning the result
  return(estimate.var)
  
}

GetConfInt <- function(list.est, alpha = 0.05) {
  # Computing the confidence interval for the estimates produced either
  # by the Ipfp() or ObtainModelEstimates() functions.
  # See confint for a new version of this functionnality.
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   list.est: The list produced by either Ipfp() or ObtainModelEstimates().
  #   alpha: The confidence level.
  #
  # Returns: a list containing the confidence interval of the estimates and
  #          their probabilities
  
  # checking that a list.est is provided
  if (is.null(list.est) == TRUE)  {
    stop('Error: a list containing the estimate and their standard deviation
         is missing!')
  }
  
  # checking that alpha is in [0,1]
  if (alpha < 0.0 | alpha > 1.0) {
    stop('Error: the confidence level alpha should be in [0,1]!')
  }
  
  # checking that list.est contains the required elements
  if (is.null(list.est$x.hat) == TRUE | is.null(list.est$x.hat.se) == TRUE) {
    stop('Error: list.est does not have the component(s) x.hat and/or x.hat.se')
  }
  if (is.null(list.est$p.hat) == TRUE | is.null(list.est$p.hat.se) == TRUE) {
    stop('Error: list.est does not have the component(s) p.hat and/or p.hat.se')
  }
  
  # checking if the lenghts of the required inputs are consistent
  if (length(list.est$x.hat) != length(list.est$x.hat.se)) {
    stop('Error: lengths of x.hat and x.hat.se components are not equal!')
  }  
  if (length(list.est$p.hat) != length(list.est$p.hat.se)) {
    stop('Error: lengths of p.hat and p.hat.se components are not equal!')
  }
  
  # computing the confidence interval
  l <- qnorm(1 - alpha * 0.5)  
  
  # ... lower bound (counts)
  ci.lo <- Array2Vector(list.est$x.hat) - l * list.est$x.hat.se
  ci.lo <- Vector2Array(ci.lo, dim(list.est$x.hat))
  dimnames(ci.lo) <- dimnames(list.est$x.hat)
  # ... upper bound (counts)
  ci.up <-Array2Vector(list.est$x.hat) + l * list.est$x.hat.se
  ci.up <- Vector2Array(ci.up, dim(list.est$x.hat))
  dimnames(ci.up) <- dimnames(list.est$x.hat)
  # ... lower bound (probabilities)  
  ci.lo.p <-Array2Vector(list.est$p.hat) - l * list.est$p.hat.se
  ci.lo.p <- Vector2Array(ci.lo.p, dim(list.est$x.hat))
  dimnames(ci.lo.p) <- dimnames(list.est$x.hat)
  # ... upper bound (probabilities)
  ci.up.p <-Array2Vector(list.est$p.hat) + l * list.est$p.hat.se
  ci.up.p <- Vector2Array(ci.up.p, dim(list.est$x.hat))
  dimnames(ci.up.p) <- dimnames(list.est$x.hat)
  
  # returning the result
  result <- list("lower.x" = ci.lo, "upper.x" = ci.up,
                 "lower.p" = ci.lo.p, "upper.p" = ci.up.p)
  return(result)
  
}
