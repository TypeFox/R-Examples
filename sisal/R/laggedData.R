### File R/laggedData.R
### This file is part of the sisal package for R.
###
### Copyright (C) 2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### A copy of the GNU General Public License is available at
### http://www.r-project.org/Licenses/

laggedData <- function(x, lags=0:9, stepsAhead=1) {
    stopifnot(is.numeric(lags), length(lags) > 0, round(lags) == lags,
              lags >= 0, is.atomic(x), !is.null(x),
              is.numeric(stepsAhead), length(stepsAhead) == 1,
              round(stepsAhead) == stepsAhead, stepsAhead >= 0,
              max(lags) + stepsAhead + 1 <= length(x))

    x2 <- as.vector(x) # drop attributes
    max.lag <- max(lags)
    n <- length(x2)
    y <- x2[(1 + max.lag + stepsAhead):n]
    n.out <- length(y)
    n.lags <- length(lags)
    X <- matrix(x2[1], n.out, n.lags)
    for(k in 1:n.lags) {
        X[, k] <- x2[seq(from = 1 + max.lag - lags[k],
                         length.out = n.out, by = 1)]
    }
    colnames(X) <- paste0("lag.", lags)
    list(X=X, y=y)
}
