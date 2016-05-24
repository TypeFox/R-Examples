controlBMAgamma <-
function(maxIter = Inf, tol = sqrt(.Machine$double.eps), 
         power = 1, startupSpeed = NULL,
         init = list(varCoefs = NULL, weights = NULL),
         optim.control = list(ndeps = rep(sqrt(.Machine$double.eps),2))) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#

 if (is.infinite(maxIter) && maxIter > 0) maxIter <- .Machine$integer.max
 list(maxIter = maxIter, tol = tol, power = power, startupSpeed = startupSpeed,
      init = init, optim.control = optim.control) 
}

