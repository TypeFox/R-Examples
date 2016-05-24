controlBMAnormal <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), equalVariance = TRUE, 
         biasCorrection = c("regression", "additive", "none"),
         init = list(sd = NULL, weights = NULL))
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, equalVariance = equalVariance,
      biasCorrection = biasCorrection[1], init = init)
}

