#  Copyright (C) 2005, 2006, 2008, 2010, 2013 David Firth and Heather Turner
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

getContrasts <- function(model, set = NULL,
                         ref = "first",
                         scaleRef = "mean",
                         scaleWeights = NULL,
                         dispersion = NULL,
                         check = TRUE,
                         ...){
    coefs <- parameters(model)
    l <- length(coefs)
    if (!l)
        stop("Model has no non-eliminated parameters")
    of.interest <- ofInterest(model)
    if (!length(of.interest)) of.interest <- seq(l)
    coefNames <- names(coefs)
    if (is.null(set))
        set <- unlist(pickFrom(coefNames[of.interest], 1, ...))
    setLength <- length(set)
    if (setLength == 0) stop(
            "No non-empty parameter set specified")
    if (setLength < 1.5) stop(
            "For contrasts, at least 2 parameters are needed in a set")
    if (is.numeric(set)) set <- coefNames[set]
    for (refName in c("ref", "scaleRef"[!is.null(scaleWeights)])) {
        refSpec <- c(get(refName))
        if (is.numeric(refSpec)){
            assign(refName, refSpec)
            if (length(refSpec) == 1){
                if (refSpec %in% seq(setLength)) {
                    temp <- rep(0, setLength)
                    temp[refSpec] <- 1
                    assign(refName, temp)
                } else stop("The specified ", refName, " is out of range")
            }
            if (length(refSpec) != setLength)
                stop("The specified ", refName, " has the wrong length")
            if ((sum(refSpec) - 1) ^ 2 > 1e-10)
                stop("The ", refName, " weights do not sum to 1")
        }
        else
            assign(refName,
                   switch(refSpec, "first" = c(1, rep.int(0, setLength - 1)),
                          "last" = c(rep.int(0, setLength - 1), 1),
                          "mean"= rep.int(1/setLength, setLength),
                          stop("Specified ", refName, " is not an opton.")))
    }
    setCoefs <- coefs[coefNames %in% set]
    contr <- setCoefs - ref %*% setCoefs
    grad <- diag(rep(1, setLength))
    grad <- grad - ref
    rownames(grad) <- set

    if (!is.null(scaleWeights)) {
        if (is.numeric(scaleWeights)) {
            scaleWeights <- c(scaleWeights)
            if (length(scaleWeights) != setLength)
                stop("The specified scaleWeights has the wrong length")
        }
        else scaleWeights <-
            switch(scaleWeights,
                   unit = rep.int(1, setLength),
                   setLength = rep.int(1/setLength, setLength),
                   stop("Specified scaleWeights is not an opton."))
        d <- setCoefs - scaleRef %*% setCoefs
        vd <- scaleWeights * d
        vdd <- sqrt(drop(vd %*% d))
        contr <- contr/vdd
        grad <- ((scaleRef * sum(vd) - vd) %o% contr/vdd + grad)/vdd
    }

    combMatrix <- matrix(0, l, setLength)
    combMatrix[match(set, coefNames), ] <- grad
    colnames(combMatrix) <- set

    Vcov <-  vcov(model, dispersion = dispersion)

    if (!is.logical(check) && !(all(check %in% seq(setLength)))) {
        stop("check must be TRUE or FALSE or a suitable numeric index vector")
    }
    iden <- rep(TRUE, ncol(combMatrix))  ## all unchecked as yet
    names(iden) <- colnames(combMatrix)
    if (is.logical(check)) {
        if (check) iden <- checkEstimable(model, combMatrix)
    }
    else iden[check] <- checkEstimable(model, combMatrix[, check])

    if (any(!na.omit(iden))) {
        if (all(!na.omit(iden))) {
            warning("None of the specified contrasts is estimable",
                    call. = FALSE)
            return(NULL)
        }
        cat("Note: the following contrasts are unestimable:\n")
        print(names(iden)[iden %in% FALSE])
    }
    not.unestimable <- iden | is.na(iden)

    combMatrix <- combMatrix[, not.unestimable, drop = FALSE]
    V <- crossprod(combMatrix, crossprod(Vcov, combMatrix))
    result <- data.frame(contr[not.unestimable], sqrt(diag(V)))
    dimnames(result) <- list(set[not.unestimable], c("Estimate", "Std. Error"))

    relerrs <- NULL
    if (sum(not.unestimable) > 2 && is.null(scaleWeights)) {
        estimable.names <- names(not.unestimable)[not.unestimable]
        Vcov <- Vcov[estimable.names, estimable.names, drop = FALSE]
        QVs <- try(qvcalc(Vcov), silent = TRUE)
        if (inherits(QVs, "try-error"))
            message("Quasi-variances could not be computed")
        else {
            quasiSE <- sqrt(QVs$qvframe$quasiVar)
            result <- cbind(result, quasiSE)
            names(result)[1:2] <- c("estimate", "SE")
            result$quasiVar <- QVs$qvframe$quasiVar
            relerrs <- QVs$relerrs
        }
    }
    return(structure(list(covmat = Vcov,
                          qvframe = result,
                          relerrs = relerrs,
                          modelcall = model$call),
                     class = "qv")
           )
}
