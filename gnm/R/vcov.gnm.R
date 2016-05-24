#  Code to estimate dispersion from summary.glm from the stats package for R.
#
#  Copyright (C) 1995-2005 The R Core Team
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

## returns vcov for the non-eliminated parameters
vcov.gnm <-  function(object, dispersion = NULL, with.eliminate = FALSE, ...){
    if (is.null(dispersion)) {
        if (any(object$family$family == c("poisson", "binomial")))
            dispersion <- 1
        else if (object$df.residual > 0) {
            if (any(object$weights == 0))
                warning("observations with zero weight ",
                        "not used for calculating dispersion")
            dispersion <- sum(object$weights * object$residuals^2)/
                object$df.residual
        }
        else dispersion <- Inf
    }
    constrain <- object$constrain
    eliminate <- object$eliminate
    nelim <- nlevels(eliminate)
    w <- as.vector(object$weights)
    X <- model.matrix(object)
    ind <- !(seq_len(ncol(X)) %in% constrain)
    cov.unscaled <- array(0, dim = rep(ncol(X), 2),
                          dimnames = list(colnames(X), colnames(X)))
    if (!length(ind)) {
        if (nelim && with.eliminate) {
            Ti <- 1/sapply(split(w, eliminate), sum)
            attr(cov.unscaled, "varElim") <- dispersion * Ti
        }
        return(structure(cov.unscaled, dispersion = dispersion,
                         ofInterest = NULL, class = "vcov.gnm"))
    }

    if (length(constrain)) X <- X[, -constrain, drop = FALSE]
    W.X <- sqrt(w) * X
    if (object$rank == ncol(W.X)) {
        cov.unscaled[ind, ind] <- chol2inv(chol(crossprod(W.X)))
    } else {
        if (is.null(eliminate)) {
            cov.unscaled[ind, ind] <- MPinv(crossprod(W.X), method = "chol",
                                            rank = object$rank)
        } else {
            ## try without ridge and generalized inverse of Q
            Ti <- 1/sapply(split(w, eliminate), sum)
            U <- rowsum(sqrt(w) * W.X, eliminate)
            W <- crossprod(W.X)
            Ti.U <- Ti * U
            UTU <- crossprod(U, Ti.U)
            cov.unscaled[ind, ind] <- MPinv(W - UTU, method = "chol",
                                            rank = object$rank - nelim)
            if (with.eliminate) {
                rownames(Ti.U) <- names(attr(coef(object), "eliminated"))
                attr(cov.unscaled, "covElim") <- dispersion *
                    -Ti.U %*% cov.unscaled[ind, ind]
                attr(cov.unscaled, "varElim") <- dispersion *
                    -rowSums(attr(cov.unscaled, "covElim") *  Ti.U) + Ti
            }
        }
    }
    structure(dispersion * cov.unscaled, dispersion = dispersion,
              ofInterest = ofInterest(object), class = "vcov.gnm")
}
