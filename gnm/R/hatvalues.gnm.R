#  Copyright (C) 2005, 2006, 2008, 2010, 2012 Heather Turner
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

hatvalues.gnm <- function(model, ...) {
    X <- as(model.matrix(model), "sparseMatrix")
    var <- unclass(vcov(model, with.eliminate = TRUE))
    eliminate <- model$eliminate
    scale <- model$weights/attr(var, "dispersion")
    hat <- rowSums((X %*% var) * X) * scale
    if (!is.null(eliminate)) {
        ## no covElim!
        if (length(model$constrain))
            X <- X[, -model$constrain, drop = FALSE]
        hat <- hat +
            (2 * rowSums(X * attr(var, "covElim")[eliminate, , drop = FALSE]) +
                 attr(var, "varElim")[eliminate]) * scale
    }
    hat <- naresid(model$na.action, hat)
    hat[is.na(hat)] <- 0
    hat[hat > 1 - 100 * .Machine$double.eps] <- 1
    if (!is.null(model$table.attr))
        attributes(hat) <- model$table.attr
    hat
}

