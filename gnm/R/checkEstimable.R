#  Copyright (C) 2005, 2006, 2008, 2010 David Firth and Heather Turner
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

checkEstimable <- function(model,
                           combMatrix = diag(length(coef(model))),
                           tolerance = NULL)
{
    if (!inherits(model, "gnm")) stop("model not of class gnm")
    coefs <- coef(model)
    l <- length(coefs)
    combMatrix <- as.matrix(combMatrix)
    if (nrow(combMatrix) != l) stop(
          "dimensions of combMatrix do not match coef(model)")
    X <- model.matrix(model)[, !is.na(coefs), drop = FALSE]
    combMatrix <- scale(combMatrix[!is.na(coefs), ], center = FALSE)
    resultNA <- apply(combMatrix, 2, function(col) any(is.na(col)))
    result <- logical(ncol(combMatrix))
    is.na(result) <- resultNA
    eliminate <- model$eliminate
    if (!is.null(eliminate)) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum(X, eliminate)/tabulate(eliminate)
        if (attr(terms(model), "intercept") == 1) subtracted[,1] <- 0
        X <- X - subtracted[eliminate, , drop = FALSE]
    }
    rankX <- model$rank - nlevels(eliminate)
    check.1 <- function(comb){
        Xc <- rbind(X, comb)
        rankXc <- quickRank(Xc, tol = tolerance)
        return(rankXc == rankX)
    }
    result[!resultNA] <- apply(combMatrix[, !resultNA, drop = FALSE],
                               2, check.1)
    names(result) <- colnames(combMatrix)
    return(result)
}
