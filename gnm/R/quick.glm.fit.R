#  Copyright (C) 2006 David Firth
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

"quick.glm.fit" <-
##  A wrapper for glm.fit, which is much faster when a large number
##  of parameters can be eliminated, but which typically (if nIter is small)
##  stops before convergence.  Useful for getting gnm starting values.
##
##  The eliminate argument is assumed numeric (no. of columns in X).  No
##  check is done here on "eliminability" of the specified columns.
##
##  The non-eliminated columns are assumed not to include the intercept (ie
##  no column of ones).
##
##  When eliminate is used, only the "coefficients" component is returned.
##  (for reasons of speed/laziness).  This is fine for gnm purposes, but if
##  quick.glm.fit is made into a `method' for glm() fits then the result
##  needs to have various other components added.
##
##  No account is taken of NAs -- will that be a problem, or have they gone by
##  the time glm.fit gets called?
##
    function (x, y,
              weights = rep(1, length(y)),
              offset = rep(0, length(y)),
              family = gaussian(),
              eliminate = 0,
              nIter = 2,
              verbose = FALSE)
{
    if (eliminate == 0)
        return(suppressWarnings(glm.fit(x, y, weights = weights,
                                        offset = offset,
                                        family = family)$coef))
    ##  The rest handles the case of eliminated columns in X
    xElim <- x[ , seq(eliminate), drop = FALSE]
    if (eliminate < ncol(x))
        xNotElim <- cbind(1, x[ , (eliminate + 1):ncol(x), drop = FALSE])
    else
        xNotElim <- matrix(1, nrow(x), 1)
    os.by.level <- numeric(eliminate)
    model <- suppressWarnings(glm.fit(xNotElim, y,
                                      weights = weights,
                                      offset = offset,
                                      family = family,
                                      control = glm.control(maxit = 1)))
    for (i in 1:nIter) {
        if (verbose) cat("quick.glm.fit iteration", i,
                         "deviance =", deviance(model), "\n")
        w <- xElim * model$weights
        wz <- w * model$residuals
        os.by.level <- os.by.level + colSums(wz)/colSums(w) + coef(model)[1]
        os.vec <- offset + colSums(os.by.level * t(xElim))
        model <- suppressWarnings(glm.fit(xNotElim, y,
                                          weights = weights,
                                          offset = os.vec,
                                          etastart = model$linear.predictors,
                                          family = family,
                                          control = glm.control(maxit = 2)))
    }
    structure(c(os.by.level + coef(model)[1], coef(model)[-1]),
              names = colnames(x))
}
