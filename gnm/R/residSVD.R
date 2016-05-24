#  Copyright (C) 2005, 2012 David Firth and Heather Turner
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

residSVD <- function(model, fac1, fac2, d = 1) {
    if (!is.null(model$call$data)) {
        Data <- as.data.frame(eval(model$call$data, parent.frame()))
        fac1 <- eval(match.call()$fac1, Data)
        fac2 <- eval(match.call()$fac2, Data)
    }
    if (!inherits(model, "glm") && class(model) != "lm") stop(
                                        "model not of class lm, glm or gnm")
    if (!is.factor(fac1)) stop("fac1 must be a factor")
    if (!is.factor(fac2)) stop("fac2 must be a factor")
    Data <- data.frame(fac1, fac2)
    if (!is.null(model$na.action)) Data <- Data[-model$na.action, ]
    weights <- if (!is.null(model$weights)) as.vector(model$weights) else 1
    X <- data.frame(rw = as.vector(model$residuals) * weights, w = weights)
    X <- lapply(X, tapply, Data, sum, simplify = TRUE)
    X <- X$rw/X$w
    X <- svd(naToZero(X), d, d)
    uPart <- sqrt(X$d[seq(d)]) * t(X$u)
    vPart <- sqrt(X$d[seq(d)]) * t(X$v)
#    uPartNegative <- apply(uPart, 1, function(row) all(row < 0))
#    vPartNegative <- apply(vPart, 1, function(row) all(row < 0))
#    multiplier <- ifelse(uPartNegative + vPartNegative == 1, -1, 1)
    multiplier <-  1
    result <- t(cbind(uPart, vPart) * multiplier)
    rownames(result) <- c(paste("fac1", levels(fac1), sep = "."),
                          paste("fac2", levels(fac2), sep = "."))
    colnames(result) <- 1:d
    drop(result)
}
