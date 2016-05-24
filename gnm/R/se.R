#  Copyright (C) 2005, 2006, 2010 David Firth and Heather Turner
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

## now only computes se for non-eliminated parameters
se <- function(model, estimate = NULL, checkEstimability = TRUE,
               Vcov = NULL, dispersion = NULL, ...){
    if (!inherits(model, "gnm")) stop("model is not of class \"gnm\"")
    if (!is.null(Vcov) && !is.null(dispersion)){
        Vcov <- Vcov * dispersion
    } else {
        Vcov <- vcov(model, dispersion = dispersion, use.eliminate = FALSE)
    }
    if (!length(Vcov)) return("Model has no non-eliminated parameters")
    coefs <- coef(model)
    coefNames <- names(coefs)
    eliminate <- model$eliminate
    nelim <- nlevels(eliminate)
    l <- length(coefs)
    if (identical(estimate, "[?]"))
        estimate <- pickCoef(model,
                             title = paste("Estimate standard errors",
                             "for one or more gnm coefficients"))
    if (is.null(estimate)){
        if (!is.null(model$ofInterest)) estimate <- ofInterest(model)
        else estimate <- seq(model$coefficients)
    }
    if (is.character(estimate))
        estimate <- match(estimate, coefNames, 0)
    if (is.vector(estimate) && all(estimate %in% seq(coefs))) {
        if (!length(estimate))
            stop("no non-eliminated coefficients specified by 'estimate'",
                 "argument")
        comb <- naToZero(coefs[estimate])
        var <- Vcov[estimate, estimate]
        coefMatrix <- matrix(0, l, length(comb))
        coefMatrix[cbind(estimate, seq(length(comb)))] <- 1
        colnames(coefMatrix) <- names(comb)
    }
    else {
        coefMatrix <- as.matrix(estimate)
        if (!is.numeric(coefMatrix))
            stop("'estimate' should specify parameters using ",
                 "\"pick\" or a vector of \n names/indices; ",
                 "or specify linear combinations using ",
                 "a numeric vector/matrix.")
        if (nrow(coefMatrix) != l)
            stop("NROW(estimate) should equal ",
                 "length(coef(model)) - nlevels(model$eliminate)")
        comb <- drop(crossprod(coefMatrix, naToZero(coefs)))
        var <- crossprod(coefMatrix, crossprod(Vcov, coefMatrix))
    }
    estimable <- rep(TRUE, ncol(coefMatrix))
    if (checkEstimability) {
        estimable <- checkEstimable(model, coefMatrix, ...)
        if (any(!na.omit(estimable)))
            cat("Std. Error is NA where estimate is fixed or unidentified\n")
    }
    if (is.matrix(var))
        sterr <- sqrt(diag(var))
    else
        sterr <- sqrt(var)
    is.na(sterr[estimable %in% c(FALSE, NA)]) <- TRUE
    result <- data.frame(comb, sterr)
    rowNames <- colnames(coefMatrix)
    if (is.null(rowNames))
        rowNames <- paste("Combination", ncol(coefMatrix))
    dimnames(result) <- list(rowNames, c("Estimate", "Std. Error"))
    result
}
