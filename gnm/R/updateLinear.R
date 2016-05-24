#  Copyright (C) 2005, 2006, 2010 Heather Turner
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

updateLinear <- function(which, theta, y, mu, eta, offset, weights, family,
                         modelTools, X, eliminate) {
    dmu <- family$mu.eta(eta)
    vmu <- family$variance(mu)
    w <- weights * dmu * dmu / vmu
    theta[which] <- 0
    offsetVarPredictors <- modelTools$varPredictors(theta)
    offset <- offset + modelTools$predictor(offsetVarPredictors)
    z <- eta - offset + (y - mu)/dmu
    if (is.null(eliminate))
        naToZero(lm.wfit(X[,which, drop = FALSE], z, w)$coef)
    else
        suppressWarnings(glm.fit.e(X[,which, drop = FALSE], z,
                                   weights = w, intercept = FALSE,
                                   eliminate = eliminate, coefonly = TRUE))
}

