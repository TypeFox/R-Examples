### File R/bootMSE.R
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

bootMSE <- function(object, dataset = NULL, R = 1000,
                    inputs = c("L.f", "L.v", "full"),
                    method = c("OLS", "magic"), standardize="inherit",
                    stepsAhead = NULL, noiseSd = NULL, verbose = 1, ...) {

    if (missing(object) || !inherits(object, "sisal")) {
        stop("'object' must exist and be of class \"sisal\"")
    }
    stopifnot(is.numeric(R), length(R) == 1, is.finite(R), R > 1,
              round(R) == R)
    if (!identical(standardize, "inherit") &&
        !identical(standardize, TRUE) && !identical(standardize, FALSE)) {
        stop("'standardize' must be \"inherit\", TRUE or FALSE")
    }
    orig.standardize <- object[["params"]][["standardize"]]
    if (is.character(standardize)) {
        standardize2 <- orig.standardize
    } else {
        standardize2 <- standardize
    }
    inputs2 <- match.arg(inputs)
    method2 <- match.arg(method)
    if (is.null(dataset)) {
        dataString <- attr(object, "dataset", exact = TRUE)
        dataDigest <- attr(object, "digest", exact = TRUE)
        if (is.null(dataString) || !is.character(dataString) ||
            length(dataString) != 1) {
            stop("must specify 'dataset'")
        } else if (grepl("laser", dataString) &&
                   identical(dataDigest, .MD5laser)) {
            dataset2 <- "laser"
        } else if (grepl("poland", dataString) &&
                   identical(dataDigest, .MD5poland.learn)) {
            dataset2 <- "poland"
        } else if (dataString %in% c("toy", "tsToy")) {
            dataset2 <- dataString
        } else {
            stop("must specify 'dataset'")
        }
    } else if (is.numeric(dataset)) {
        theDim <- dim(dataset)
        if (!is.null(theDim) && (length(theDim) > 2 || sum(theDim > 1) > 1)) {
            stop("numeric 'dataset' must be a vector or a matrix with 1 row or column")
        }
        dataset2 <- "ts"
        tsData <- as.vector(dataset)
    } else if(is.list(dataset)) {
        if (!all(c("X", "y") %in% names(dataset))) {
            stop("'dataset' of type list must have elements \"X\" and \"y\"")
        }
        X <- dataset[["X"]]
        y <- dataset[["y"]]
        dataset2 <- "X, y"
    } else {
        dataset2 <- match.arg(dataset, c("laser", "poland", "toy", "tsToy"))
    }
    if (dataset2 == "ts") {
        foo <- digestVector(dataset)
        if (identical(foo, .MD5laser.cont) ||
            (is.double(dataset) &&
             identical(round(dataset), dataset) &&
             identical(digestVector(as.integer(dataset)), .MD5laser.cont))) {
            dataset2 <- paste(dataset2, "(laser.cont)")
        } else if (identical(foo, .MD5poland.test)) {
            dataset2 <- paste(dataset2, '(poland[["test"]])')
        }
    }
    if (verbose > 0) {
        splitcat(gettextf('Dataset is "%s"', dataset2, domain="R-sisal"))
    }
    mean.X <- object[["mean.X"]]
    mean.y <- object[["mean.y"]]
    sd.X <- object[["sd.X"]]
    sd.y <- object[["sd.y"]]
    divide.y <- !object[["zeroRange.y"]] && sd.y != 0
    constant.X <- object[["constant.X"]]
    doSubset <- inputs2 %in% c("L.f", "L.v")
    if (doSubset) {
        subInputs <- object[[inputs2]]
        mean.X <- mean.X[subInputs]
        sd.X <- sd.X[subInputs]
        constant.X <- constant.X[subInputs]
    }
    if (dataset2 == "toy") {
        if (is.null(noiseSd)) {
            if (verbose >= 2) {
                splitcat(gettextf("Copying value of '%s' from 'object'",
                                  "noiseSd", domain="R-sisal"))
            }
            noiseSd2 <- attr(object, "noiseSd", exact = TRUE)
        } else {
            stopifnot(is.numeric(noiseSd), length(noiseSd) == 1, noiseSd >= 0)
            noiseSd2 <- noiseSd
        }
        if (verbose > 0) {
            splitcat(gettextf("See ?%s for a description of the data",
                              "toy.test", domain="R-sisal"))
            if (noiseSd2 == 0) {
                splitcat(gettext("Not adding noise", domain="R-sisal"))
            } else {
                splitcat(gettextf("Adding noise with standard deviation %f",
                                  noiseSd2, domain="R-sisal"))
            }
        }
        toy.test <- sisal::toy.test
        cnames <- colnames(toy.test)
        X <- toy.test[, grep("^X", cnames), drop=FALSE]
        y <- toy.test[, cnames == "y"] +
            noiseSd2 * toy.test[, cnames == "noise"]
        if (doSubset) {
            X <- X[, object[[inputs2]], drop=FALSE]
        }
    } else if (dataset2 != "X, y") {
        if (verbose > 0 &&
            (dataset2 == "laser" || dataset2 == "poland")) {
            splitcat(gettext("Downloading data", domain="R-sisal"))
        }
        if (is.null(stepsAhead)) {
            if (verbose >= 2) {
                splitcat(gettextf("Copying value of '%s' from 'object'",
                                  "stepsAhead", domain="R-sisal"))
            }
            stepsAhead2 <- attr(object, "stepsAhead", exact = TRUE)
        } else {
            stopifnot(is.numeric(stepsAhead), length(stepsAhead) == 1)
            stepsAhead2 <- stepsAhead
        }
        if (dataset2 == "laser") {
            tsData <- sisalData("laser.cont")
            if (verbose > 0) {
                cat("laser.cont\n")
            }
        } else if (dataset2 == "poland") {
            tsData <- sisalData("poland")[["test"]]
            if (verbose > 0) {
                cat('poland[["test"]]\n')
            }
        }  else if (dataset2 == "tsToy") {
            tsData <- sisal::tsToy.test
            if (verbose > 0) {
                cat("tsToy.test\n")
            }
        }
        lags <- as.integer(sub("lag\\.(\\d+)", "\\1",
                               object[["var.names"]]))
        if (doSubset) {
            lags <- lags[object[[inputs2]], drop=FALSE]
        }
        if (verbose > 0) {
            if (grepl("laser|poland", dataset2)) {
                splitcat(gettextf("See ?%s for a description of the data",
                                  "sisalData", domain="R-sisal"))
            }
            splitcat(sprintf(ngettext(stepsAhead2,
                                      "Predicting %.0f step ahead with lags:",
                                      "Predicting %.0f steps ahead with lags:",
                                      domain="R-sisal"),
                             stepsAhead2))
            print(lags)
        }
        if (length(lags) > 0) {
            ld <- laggedData(x = tsData, lags = lags, stepsAhead = stepsAhead2)
            X <- ld[["X"]]
            y <- ld[["y"]]
        } else {
            y <- tsData
            X <- matrix(numeric(0), length(y), 0)
        }
    }
    if (length(y) != nrow(X)) {
        stop("mismatch in sizes of 'X' and 'y'")
    }
    if (standardize2) {
        y <- y - mean.y
        if (divide.y) {
            y <- y / sd.y
        }
    }

    coefs <- NULL
    if (method2 == "OLS") {
        theModel <- object[[paste0("lm.", inputs2)]]
        if (!is.null(theModel)) {
            coefs <- coef(theModel)
        }
    } else {
        theModel <- object[[paste0("magic.", inputs2)]]
        if (!is.null(theModel)) {
            coefs <- theModel[["b"]]
        }
    }
    if (!is.null(coefs)) {
        coefNames <- names(coefs)
        if (is.null(coefNames)) {
            stop("missing names")
        }
        iceptFlag <- coefNames != "(Intercept)"
        coefNames <- coefNames[iceptFlag]
        coefs <- coefs[iceptFlag]
        colNames <- colnames(X)
        if (is.null(colNames) || !all(coefNames %in% colNames)) {
            stop("missing variables or names")
        }
        if (anyDuplicated(coefNames) || anyDuplicated(colNames)) {
            stop("duplicated names")
        }
        X <- X[, match(coefNames, colNames), drop = FALSE]
        NA.flag <- rowSums(is.na(X)) > 0 | is.na(y)
        n.NA <- sum(NA.flag)
        if (n.NA > 0) {
            warning(gettextf("dropping %d missing samples out of a total of %d",
                             n.NA, length(y), domain = "R-sisal"), domain = NA)
            notNA <- which(!NA.flag)
            X <- X[notNA, , drop = FALSE]
            y <- y[notNA]
        }
        if (length(y) == 0) {
            stop("no samples available")
        }
        mean.X.new <- colMeans(X)
        zeroRange.X.new <- logical(ncol(X))
        for (k in seq_along(zeroRange.X.new)) {
            zeroRange.X.new[k] <- zeroRange(X[, k], mean.x = mean.X.new[k])
        }
        zero.X.new <- zeroRange.X.new
        zero.X.new[zero.X.new] <- vapply(mean.X.new[zero.X.new], function (x) {
            isTRUE(all.equal(0, x, check.attributes = FALSE))
        }, TRUE)
        constant.X.new <- zeroRange.X.new & !zero.X.new
        same.mean <- vapply(mean.X.new - mean.X,
                            function(x) isTRUE(all.equal(0, x)),
                            logical(1))
        not.constant.all <- !(constant.X & constant.X.new & same.mean)
        X <- X[, not.constant.all, drop=FALSE]
        coefs <- coefs[not.constant.all]
        sd.X <- sd.X[not.constant.all]
        mean.X <- mean.X[not.constant.all]
        constant.X <- constant.X[not.constant.all]
        sd.X.tmp <- sd.X
        zeroRange.X <- object[["zeroRange.X"]]
        if (doSubset) {
            zeroRange.X <- zeroRange.X[subInputs]
        }
        zeroRange.X <- zeroRange.X[not.constant.all]
        sd.X.tmp[zeroRange.X | sd.X.tmp == 0] <- 1
        if (standardize2) {
            mean.X.tmp <- mean.X
            mean.X.tmp[constant.X] <- mean.X[constant.X] - 1
            X <- sweep(sweep(X, 2, mean.X.tmp, "-"), 2, sd.X.tmp, "/")
            if (!orig.standardize) {
                coefs <- coefs * sd.X.tmp
                if (divide.y) {
                    coefs <- coefs / sd.y
                }
            }
        } else {
            X <- cbind(1, X)
            if (orig.standardize) {
                coefs <- coefs / sd.X.tmp
                if (divide.y) {
                    coefs <- coefs * sd.y
                }
            }
            icept <- mean.y - drop(crossprod(coefs, mean.X))
            coefs <- c("(Intercept)"=icept, coefs)
        }
        predicted <- drop(X %*% coefs)
        error <- y - predicted
    } else {
        warning("no variables in the model")
        NA.flag <- is.na(y)
        n.NA <- sum(NA.flag)
        if (n.NA > 0) {
            warning(gettextf("dropping %d missing samples out of a total of %d",
                             n.NA, length(y), domain = "R-sisal"), domain = NA)
            error <- y[!NA.flag]
        } else {
            error <- y
        }
        if (length(error) == 0) {
            stop("no samples available")
        }
    }

    cl <- as.call(c(as.name("boot"), list(...)))
    fArgs <- lapply(as.list(match.call(boot, cl))[-1L],
                     function(x) call("quote", x))
    fArgs[c("data", "statistic")] <- NULL
    argNames <- names(fArgs)
    fArgs[["R"]] <- R
    fArgs[["stype"]] <- "i"
    if (!"sim" %in% argNames) {
        fArgs[["sim"]] <- "ordinary"
    }
    if (!"m" %in% argNames) {
        fArgs[["m"]] <- 0
    }
    if (!"simple" %in% argNames && fArgs[["sim"]] == "ordinary" &&
        is.numeric(fArgs[["m"]]) && sum(fArgs[["m"]]) == 0) {
        fArgs[["simple"]] <- TRUE
    }
    do.call("boot", c(fArgs, alist(data = error^2,
                                   statistic = function(data, idx) {
                                       mean(data[idx])
                                   })))
}
