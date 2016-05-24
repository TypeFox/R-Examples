#  Copyright (C) 2005, 2008, 2010, 2012, 2014, 2015 Heather Turner
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

predict.gnm <- function (object, newdata = NULL,
                         type = c("link", "response", "terms"), se.fit = FALSE,
                         dispersion = NULL, terms = NULL,
                         na.action = na.exclude, ...) {
    type <- match.arg(type)
    if (type == "terms") {
        hasintercept <- attr(object$terms, "intercept") > 0L
        ## do not include eliminate term - cannot check estimability without
        ## creating full matrix, defeating point of eliminate
        if (is.null(terms)) {
            terms <- attr(object$terms, "term.labels")
        } else {
            terms <- setdiff(terms, "(eliminate)")
        }
    }
    if (missing(newdata)) {
        pred <- switch(type, link = object$predictors,
                       response = object$fitted.values,
                       terms = {pred <- termPredictors(object)
                                ## see 6.3.6 white book & predict.lm
                                if (hasintercept) {
                                    predc <- sweep(pred, 2, colMeans(pred))
                                    const <- sum(pred[1,]) - sum(predc[1,])
                                    structure(predc[, terms, drop = FALSE],
                                              constant = const)
                                } else structure(pred[, terms, drop = FALSE],
                                                 constant = 0)})
        if (!is.null(na.act <- object$na.action)){
            pred <- napredict(na.act, pred)
        }
        if (!inherits(pred, "matrix") && !is.null(object$table.attr))
            attributes(pred) <- object$table.attr
    } else {
        modelTerms <- delete.response(terms(object))
        ## evaluate eliminate in environment of formula
        if (is.null(object$eliminate)){
            modelData <- model.frame(modelTerms, newdata, na.action = na.action,
                                     xlev = object$xlevels)
        } else {
            modelData <-
                model.frame(modelTerms, newdata, eliminate = eval(eliminate),
                            na.action = na.action, xlev = object$xlevels)
        }
        ## use same contrasts as in original model
        contr <- lapply(model.frame(object)[names(modelData)],
                        attr, "contrasts")
        for (i in which(!sapply(contr, is.null))){
            modelData[[i]] <- C(modelData[[i]], contr[[i]])
        }
        if (length(offID <- attr(modelTerms, "offset"))){
            offset <- eval(attr(modelTerms, "variables")[[offID + 1]],
                           newdata)
        } else offset <- eval(object$call$offset, newdata)
        modelTools <- gnmTools(modelTerms, modelData)
        varPredictors <- modelTools$varPredictors(parameters(object))
        pred <- modelTools$predictor(varPredictors, term = type == "terms")
        if (type == "terms") {
            rownames(pred) <- rownames(modelData)
        } else names(pred) <- rownames(modelData)
        if (!is.null(offset))  pred <- offset + pred
        if (!is.null(object$eliminate)) {
            prede <- attr(coef(object), "eliminate")
            if (type != "terms") pred <- prede[modelData$`(eliminate)`] + pred
        }
        switch(type, response = {pred <- family(object)$linkinv(pred)},
               terms = {if (hasintercept) {
                            predc <- sweep(pred, 2,
                                           colMeans(termPredictors(object)))
                            const <- sum(pred[1,]) - sum(predc[1,])
                            pred <- structure(predc[, terms, drop = FALSE],
                                              constant = const)
                        } else structure(pred[, terms, drop = FALSE],
                                         constant = 0)},
               link = )
        if (!is.null(na.act <- attr(modelData, "na.action")))
            pred <- napredict(na.act, pred)
    }
    if (se.fit) {
        V <- vcov(object, dispersion = dispersion, with.eliminate = TRUE)
        residual.scale <- as.vector(sqrt(attr(V, "dispersion")))
        if (missing(newdata)) {
            X <- model.matrix(object)
            elim <- object$eliminate
        } else {
            X <- modelTools$localDesignFunction(parameters(object),
                                                varPredictors)
            elim <- modelData$`(eliminate)`
        }
        covElim <- attr(V, "covElim")[elim, , drop = FALSE]
        varElim <- attr(V, "varElim")[elim]
        switch(type,
               link = {
                   if (is.null(elim))
                       se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))
                   else se.fit <-
                       sqrt(diag(X %*% tcrossprod(V, X)) +
                                2 * rowSums(X * covElim) + varElim)},
               response = {
                   eta <- na.omit(c(family(object)$linkfun(pred)))
                   d <- family(object)$mu.eta(eta)
                   if (is.null(object$eliminate))
                       se.fit <- sqrt(diag(X %*% tcrossprod(V, X)))
                   else se.fit <- sqrt(diag(X %*% tcrossprod(V, X)) +
                                       2*rowSums(X * covElim) + varElim)
                   se.fit <- d * se.fit},
               terms = {
                   if (missing(newdata)) {
                       assign <- split(seq(ncol(X)), attr(X, "assign"))
                   } else {
                       M <- model.matrix(object)
                       assign <- split(seq(ncol(X)), attr(M, "assign"))
                   }
                   if (hasintercept) {
                       if (missing(newdata)) {
                           X <- sweep(X, 2, colMeans(X))
                       } else X <- sweep(X, 2, colMeans(M))
                   }
                   se.fit <- matrix(, nrow = nrow(X), ncol = length(terms))
                   s <- 0
                   adj <- hasintercept
                   for (i in match(terms, colnames(pred))) {
                       s <- s + 1
                       t <- assign[[i + adj]]
                       se.fit[, s] <-
                           sqrt(diag(X[, t] %*%
                                         tcrossprod(V[t, t], X[, t])))
                       ## check estimability of term
                       Xt <- X
                       Xt[, -t] <- 0
                       estimable <- checkEstimable(object, t(Xt))
                       is.na(se.fit)[estimable %in% c(FALSE, NA), s] <- TRUE
                   }
               })
        ## check estimability of predictions
        if (!missing(newdata) && type != "terms"){
            estimable <- checkEstimable(object, t(X))
            is.na(se.fit)[estimable %in% c(FALSE, NA)] <- TRUE
        }
        if (!is.null(na.act)) {
            se.fit <- napredict(na.act, se.fit)
        }
        if (inherits(pred, "table"))
            attributes(se.fit) <- object$table.attr
        else
            attributes(se.fit) <- attributes(pred)
        pred <- list(fit = pred, se.fit = se.fit,
                     residual.scale = residual.scale)
    }
    pred
}
