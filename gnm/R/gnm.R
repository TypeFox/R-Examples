#  Designed to take similar arguments to glm from the stats package from R;
#  some of the code to handle the arguments is copied/modified from glm.
#
#  Copyright (C) 1995-2005 The R Core Team
#  Copyright (C) 2005-2010, 2012, 2013 Heather Turner and David Firth
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

gnm <- function(formula, eliminate = NULL, ofInterest = NULL,
                constrain = numeric(0), #index of non-eliminated parameters
                constrainTo = numeric(length(constrain)), family = gaussian,
                data = NULL, subset, weights, na.action,  method = "gnmFit",
                checkLinear = TRUE, offset, start = NULL,
                etastart = NULL, mustart = NULL, tolerance = 1e-6,
                iterStart = 2, iterMax = 500, trace = FALSE, verbose = TRUE,
                model = TRUE, x = TRUE, termPredictors = FALSE,
                ridge = 1e-8, ...) {

    call <- match.call()

    modelTerms <- gnmTerms(formula, substitute(eliminate), data)
    modelData <- as.list(match.call(expand.dots = FALSE))
    if (inherits(data, "table") && missing(na.action))
        modelData$na.action <- "na.exclude"
    argPos <- match(c("eliminate", "data", "subset", "weights", "na.action",
                      "offset", "etastart", "mustart"),
                    names(modelData), 0)
    modelData <- as.call(c(as.name("model.frame"),
                           formula = modelTerms,
                           modelData[argPos],
                           drop.unused.levels = TRUE))
    modelData <- eval(modelData, parent.frame())

    if (!missing(eliminate)) {
        eliminate <- modelData$`(eliminate)`
        if (!is.factor(eliminate))
            stop("'eliminate' must be a factor")
        xtf <- xtfrm(modelData$`(eliminate)`)
        ord <- order(xtf)
        if (ordTRUE <- !identical(ord, xtf)) {
            modelData <- modelData[ord, , drop = FALSE]
            eliminate <- modelData$`(eliminate)`
        }
        nElim <- nlevels(eliminate)
    }
    else nElim <- 0

    if (method == "model.frame")
        return(modelData)
    else if (!method %in% c("gnmFit", "coefNames", "model.matrix") &&
             !is.function(get(method))) {
        warning("function ", method, " can not be found. Using \"gnmFit\".\n",
                call. = FALSE)
        method <- "gnmFit"
    }

    nobs <- nrow(modelData)
    y <- model.response(modelData, "numeric")
    if (is.null(y))
        y <- rep(0, nobs)

    weights <- as.vector(model.weights(modelData))
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights are not allowed")
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    offset <- as.vector(model.offset(modelData))
    if (is.null(offset))
        offset <- rep.int(0, nobs)

    mustart <- model.extract(modelData, "mustart")
    etastart <- model.extract(modelData, "etastart")

    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("`family' not recognized")
    }

    if (family$family == "binomial") {
        if (is.factor(y) && NCOL(y) == 1)
            y <- y != levels(y)[1]
        else if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
        }
    }

    if (is.empty.model(modelTerms) && missing(eliminate)) {
        if (method == "coefNames") return(numeric(0))
        else if (method == "model.matrix")
            return(model.matrix(modelTerms, data = modelData))
        if (!family$valideta(offset))
            stop("invalid predictor values in empty model")
        mu <- family$linkinv(offset)
        if (!family$validmu(mu))
            stop("invalid fitted values in empty model")
        dmu <- family$mu.eta(offset)
        dev <- sum(family$dev.resids(y, mu, weights))
        modelAIC <- suppressWarnings(family$aic(y, rep.int(1, nobs), mu,
                                                weights, dev))
        fit <- list(coefficients = numeric(0), constrain = numeric(0),
                    constrainTo = numeric(0), eliminate = NULL,
                    predictors = offset, fitted.values = mu, deviance = dev,
                    aic = modelAIC, iter = 0,
                    weights = weights*dmu^2/family$variance(mu),
                    residuals = (y - mu)/dmu, df.residual = nobs, rank = 0,
                    family = family, prior.weights = weights, y = y,
                    converged = NA)
        if (x) fit <- c(fit, x = model.matrix(modelTerms, data = modelData))
        if (termPredictors) fit <- c(fit, termPredictors =
                                     matrix(, nrow(modelData), 0))
    }
    else {
        onlyLin <- checkLinear && all(attr(modelTerms, "type") == "Linear")
        if (onlyLin) {
            if (nElim) {
                X <- model.matrix(update(modelTerms, . ~ . + 1), modelData)
                asgn <- attr(X, "assign")
                X <- X[,-1, drop = FALSE]
                attr(X, "assign") <- asgn[-1]
            }
            else X <- model.matrix(modelTerms, modelData)
            coefNames <- colnames(X)
        }
        else {
            modelTools <- gnmTools(modelTerms, modelData,
                                   method == "model.matrix" | x)
            coefNames <- names(modelTools$start)
        }
        if (method == "coefNames") return(coefNames)
        nParam <- length(coefNames)

        if (identical(constrain, "[?]"))
            call$constrain <- constrain <-
                unlist(pickFrom(coefNames,
                                edit.setlabels = FALSE,
                                title =
                                "Constrain one or more gnm coefficients",
                                items.label = "Model coefficients:",
                                warningText =
                                "No parameters were specified to constrain",
                                return.indices = TRUE))
        if (is.character(constrain)) {
            res <- match(constrain, coefNames, 0)
            if (res == 0 && length(constrain) == 1){
                constrain <- match(grep(constrain, coefNames),
                                       seq_len(nParam),  0)
            }
            else constrain <- res
        }
        ## dropped logical option
        if (!all(constrain %in% seq_len(nParam)))
            stop(" cannot match 'constrain' to non-eliminated parameters. ")

        if (is.null(start))
            start <- rep.int(NA, nElim + nParam)
        else if (length(start) != nElim + nParam) {
            if (!missing(eliminate) && length(start) == nParam)
                start <- c(rep.int(NA, nElim), start)
            else
                stop("length(start) must either equal the no. of parameters\n",
                     "or the no. of non-eliminated parameters.")
        }

        if (onlyLin) {
            if (length(constrain)) {
                offset <- drop(offset +
                               X[, constrain, drop = FALSE] %*% constrainTo)
                X[, constrain] <- 0
            }
            if (method == "model.matrix") return(X)
        }
        else if (method == "model.matrix"){
            theta <- modelTools$start
            theta[!is.na(start)] <- start[!is.na(start)]
            theta[constrain] <- constrainTo
            theta[is.na(theta)] <- seq(start)[is.na(theta)]
            varPredictors <- modelTools$varPredictors(theta)
            X <- modelTools$localDesignFunction(theta, varPredictors)
            attr(X, "assign") <- modelTools$termAssign
            return(X)
        }

        if (!is.numeric(tolerance) || tolerance <= 0)
            stop("value of 'tolerance' must be > 0")
        if (!is.numeric(iterMax) || iterMax < 0)
            stop("maximum number of iterations must be >= 0")

        if (onlyLin) {
            if (any(is.na(start))) start <- NULL
            fit <- glm.fit.e(X, y, weights = weights, start = start,
                             etastart = etastart, mustart = mustart,
                             offset = offset, family = family,
                             control = glm.control(tolerance, iterMax, trace),
                             intercept = attr(modelTerms, "intercept"),
                             eliminate = eliminate)
            if (sum(is.na(coef(fit))) > length(constrain)) {
                extra <- setdiff(which(is.na(coef(fit))), constrain)
                ind <- order(c(constrain, extra))
                constrain <- c(constrain, extra)[ind]
                constrainTo <- c(constrainTo, numeric(length(extra)))[ind]
            }
            if (!is.null(fit$null.deviance)) {
                extra <- match(c("effects",  "R", "qr", "null.deviance",
                                 "df.null", "boundary"), names(fit))
                fit <- fit[-extra]
            }
            names(fit)[match("linear.predictors", names(fit))] <- "predictors"
            fit$constrain <- constrain
            fit$constrainTo <- constrainTo
            if (x) {
                fit$x <- X
            }
            if (termPredictors) {
                modelTools <- gnmTools(modelTerms, modelData)
                varPredictors <- modelTools$varPredictors(naToZero(coef(fit)))
                fit$termPredictors <- modelTools$predictor(varPredictors,
                                                           term = TRUE)
            }
        }
        else if (method != "gnmFit")
            fit <- do.call(method, list(modelTools = modelTools, y = y,
                                        constrain = constrain,
                                        constrainTo = constrainTo,
                                        eliminate = eliminate, family = family,
                                        weights = weights, offset = offset,
                                        nobs = nobs, start = start,
                                        etastart = etastart, mustart = mustart,
                                        tolerance = tolerance,
                                        iterStart = iterStart,
                                        iterMax = iterMax, trace = trace,
                                        verbose = verbose, x = x,
                                        termPredictors = termPredictors,
                                        ridge = ridge, ...))
        else
            fit <- gnmFit(modelTools = modelTools, y = y, constrain = constrain,
                          constrainTo = constrainTo, eliminate = eliminate,
                          family = family, weights = weights, offset = offset,
                          nobs = nobs, start = start, etastart = etastart,
                          mustart = mustart, tolerance = tolerance,
                          iterStart = iterStart, iterMax = iterMax,
                          trace = trace, verbose = verbose, x = x,
                          termPredictors = termPredictors,
                          ridge = ridge)
    }

    if (is.null(fit)) {
        warning("Algorithm failed - no model could be estimated", call. = FALSE)
        return()
    }

    if (is.null(ofInterest) && !missing(eliminate))
        ofInterest <- seq_len(nParam)
    if (identical(ofInterest, "[?]"))
        call$ofInterest <- ofInterest <-
            pickCoef(fit,
                     warningText = paste("No subset of coefficients selected",
                     "- assuming all are of interest."))
    if (is.character(ofInterest)) {
        if (length(ofInterest) == 1)
            ofInterest <- match(grep(ofInterest, coefNames), seq_len(nParam), 0)
        else
            ofInterest <- match(ofInterest, coefNames, 0)
        if (!sum(ofInterest)) ofInterest <- seq_len(nParam)
    }
    if (!is.null(ofInterest)) {
        if (!all(ofInterest %in% seq_len(nParam)))
            stop("'ofInterest' does not specify a subset of the ",
                 "non.eliminated coefficients.")
        names(ofInterest) <- coefNames[ofInterest]
    }

    if (missing(data))
        data <- environment(formula)
    fit <- c(list(call = call, formula = formula,
                  terms = modelTerms, data = data, eliminate = eliminate,
                  ofInterest = ofInterest,
                  na.action = attr(modelData, "na.action"),
                  xlevels = .getXlevels(modelTerms, modelData),
                  offset = offset, tolerance = tolerance, iterStart = iterStart,
                  iterMax = iterMax), fit)

    if (!missing(eliminate) && ordTRUE) {
        reorder <- order(ord)
        fit <- within(fit, {
            y <- y[reorder]
            fitted.values <- fitted.values[reorder]
            predictors <- predictors[reorder]
            residuals <- residuals[reorder]
            weights <- weights[reorder]
            prior.weights <- prior.weights[reorder]
            eliminate <- eliminate[reorder]
            offset <- offset[reorder]
        })
        modelData <- modelData[reorder, , drop = FALSE]
        y <- y[reorder]
        if (x) {
            asgn <- attr(fit$x, "assign")
            fit$x <- fit$x[reorder, , drop = FALSE]
            attr(fit$x, "assign") <- asgn
        }
    }

    asY <- c("predictors", "fitted.values", "residuals", "prior.weights",
             "weights", "y", "offset")
    if (inherits(data, "table") &&
        (is.null(fit$na.action) | inherits(fit$na.action, "exclude"))) {
        attr <- attributes(data)
        if (!missing(subset)) {
            ind <- as.numeric(names(y))
            lev <- do.call("expand.grid", attr$dimnames)[ind,, drop = FALSE]
            attr$dimnames <- apply(lev, 2, unique)
            attr$dim <- unname(sapply(attr$dimnames, length))
        }
        fit$table.attr <- attr
    }
    fit[asY] <- lapply(fit[asY], structure, dim = NULL, names = names(y))
    if (termPredictors)
        rownames(fit$termPredictors) <- names(y)
    if (model)
        fit$model <- modelData
    class(fit) <- c("gnm", "glm", "lm")
    attr(fit, ".Environment") <- environment(gnm)
    fit
}

