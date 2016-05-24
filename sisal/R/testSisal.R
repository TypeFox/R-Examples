### File R/testSisal.R
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

testSisal <- function(dataset = c("tsToy", "laser", "poland", "toy"),
                      nData = Inf, FUN = "sisal", lags = NULL,
                      stepsAhead = 1, noiseSd = 0.2, verbose = 1, ...) {
    POLAND_LAGS <- 0:14
    LASER_LAGS <- 0:19
    TSTOY_LAGS <- 0:9
    FUN2 <- match.fun(FUN)
    if (is.numeric(dataset)) {
        theDim <- dim(dataset)
        if (!is.null(theDim) && (length(theDim) > 2 || sum(theDim > 1) > 1)) {
            stop("numeric 'dataset' must be a vector or a matrix with 1 row or column")
        }
        dataset2 <- "ts"
        tsData <- as.vector(dataset)
    } else {
        dataset2 <- match.arg(dataset)
        stopifnot(is.numeric(noiseSd), length(noiseSd) == 1, noiseSd >= 0)
    }
    stopifnot(is.numeric(nData), length(nData) == 1, nData > 1,
              round(nData) == nData)
    lags2 <- lags
    if (dataset2 == "ts") {
        foo <- digestVector(dataset)
        if (identical(foo, .MD5laser) ||
            (is.double(dataset) &&
             identical(round(dataset), dataset) &&
             identical(digestVector(as.integer(dataset)), .MD5laser))) {
            dataset2 <- paste(dataset2, "(laser)")
            if (is.null(lags)) {
                lags2 <- LASER_LAGS
            }
        } else if (identical(foo, .MD5poland.learn)) {
            dataset2 <- paste(dataset2, "(poland)")
            if (is.null(lags)) {
                lags2 <- POLAND_LAGS
            }
        } else if (is.null(lags)) {
            stop("'lags' must be given when dataset is a generic numeric vector")
        }
    }
    if (is.function(FUN)) {
        foo <- match.call()[["FUN"]]
        if (is.name(foo)) {
            funName <- as.character(foo)
        } else if (is.call(foo)) {
            char1 <- as.character(foo[[1]])
            if (char1 %in% c("::", ":::")) {
                funName <- gettextf("%s (package %s)",
                                    as.character(foo[[3]]), sQuote(foo[[2]]),
                                    domain="R-sisal")
            } else if (char1 == "function") {
                funName <- gettext("anonymous function", domain="R-sisal")
            } else {
                funName <- "??"
            }
        }
    } else if (is.name(FUN)) {
        funName <- as.character(FUN)
    } else if (is.character(FUN) && length(FUN) == 1) {
        funName <- FUN
    } else {
        funName <- "??"
    }
    if (verbose > 0) {
        splitcat(gettextf('Dataset is "%s"', dataset2, domain="R-sisal"))
    }
    if (dataset2 == "toy") {
        if (verbose > 0) {
            splitcat(gettextf("See ?%s for a description of the data",
                              "toy.learn", domain="R-sisal"))
            if (noiseSd == 0) {
                splitcat(gettext("Not adding noise", domain="R-sisal"))
            } else {
                splitcat(gettextf("Adding noise with standard deviation %f",
                                  noiseSd, domain="R-sisal"))
            }
        }
        toy.learn <- sisal::toy.learn
        cnames <- colnames(toy.learn)
        X <- toy.learn[, grep("^X", cnames), drop=FALSE]
        y <- toy.learn[, cnames == "y"] +
            noiseSd * toy.learn[, cnames == "noise"]
    } else {
        if (verbose > 0 &&
            (dataset2 == "laser" || dataset2 == "poland")) {
            splitcat(gettext("Downloading data", domain="R-sisal"))
        }
        if (dataset2 == "laser") {
            if (is.null(lags)) {
                lags2 <- LASER_LAGS
            }
            tsData <- sisalData("laser")
            if (verbose > 0) {
                cat("laser\n")
            }
        } else if (dataset2 == "poland") {
            if (is.null(lags)) {
                lags2 <- POLAND_LAGS
            }
            tsData <- sisalData("poland")[["learn"]]
            if (verbose > 0) {
                cat('poland[["learn"]]\n')
            }
        }  else if (dataset2 == "tsToy") {
            if (is.null(lags)) {
                lags2 <- TSTOY_LAGS
            }
            tsData <- sisal::tsToy.learn
            if (verbose > 0) {
                cat("tsToy.learn\n")
            }
        }
        if (verbose > 0) {
            if (grepl("laser|poland", dataset2)) {
                splitcat(gettextf("See ?%s for a description of the data",
                                  "sisalData", domain="R-sisal"))
            } else if (dataset2 == "tsToy") {
                splitcat(gettextf("See ?%s for a description of the data",
                                  "tsToy.learn", domain="R-sisal"))
            }
            splitcat(sprintf(ngettext(stepsAhead,
                                      "Predicting %.0f step ahead with lags:",
                                      "Predicting %.0f steps ahead with lags:",
                                      domain="R-sisal"),
                             stepsAhead))
            print(lags2)
        }
        ld <- laggedData(x = tsData, lags = lags2, stepsAhead = stepsAhead)
        X <- ld[["X"]]
        y <- ld[["y"]]
    }

    nMax <- length(y)
    if (nData < nMax) {
        if (verbose > 0) {
            splitcat(gettextf("Using sample size %.0f / %.0f", nData, nMax,
                              domain="R-sisal"))
        }
        idx <- sample.int(n = nMax, size = nData, replace = FALSE)
        X <- X[idx, , drop=FALSE]
        y <- y[idx]
    }

    cl <- as.call(c(FUN2, list(...)))
    fArgs <- lapply(as.list(match.call(FUN2, cl))[-1L],
                    function(x) call("quote", x))
    fArgs[["verbose"]] <- verbose
    fArgs[["X"]] <- quote(X)
    fArgs[["y"]] <- quote(y)
    if (verbose > 0) {
        splitcat(gettextf("Calling %s", funName, domain="R-sisal"))
        cat("\n")
    }
    res <- do.call(FUN2, fArgs)
    attr(res, "dataset") <- dataset2
    if (exists("tsData", inherits = FALSE)) {
        attr(res, "digest") <- digestVector(tsData)
        attr(res, "stepsAhead") <- stepsAhead
    }
    if (dataset2 == "toy") {
        attr(res, "noiseSd") <- noiseSd
    }
    res
}
