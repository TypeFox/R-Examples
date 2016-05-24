mdm <- function (formula, data, weights, subset, na.action,
    MaxNWts, maxit = 1000, contrasts = NULL, Hess = FALSE,
    censored = FALSE, model = TRUE, use.shortcut = TRUE, ...)
{

# file nnet/multinom.R
# copyright (C) 1994-2006 W. N. Venables and B. D. Ripley
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
#

# file MDM/mdm.R
# mdm is modified version of multinom from nnet package by B.R.Ripley
# and Bill Venables

    class.ind <- function(cl) {
        n <- length(cl)
        x <- matrix(0, n, length(levels(cl)))
        x[(1L:n) + n * (as.vector(unclass(cl)) - 1L)] <- 1
        dimnames(x) <- list(names(cl), levels(cl))
        x
    }
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$MaxNWts <- m$maxit <- m$summ <- m$Hess <- m$contrasts <- m$censored <- m$q <- m$model <- m$use.shortcut <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    X <- model.matrix(Terms, m, contrasts)
    cons <- attr(X, "contrasts")
    Xr <- qr(X)$rank
    Y <- model.response(m)
    w <- model.weights(m)
    if (length(w) == 0L)
        if (is.matrix(Y))
            w <- rep(1, dim(Y)[1L])
        else w <- rep(1, length(Y))
    offset <- model.offset(m)
    r <- ncol(X)
    if (is.matrix(Y)) {
		use.shortcut <-  use.shortcut & ((Xr==1)|(Xr==nrow(X)))
		p <- ncol(Y)
        sY <- Y %*% rep(1, p)
        if (any(sY == 0))
            stop("some case has no observations")
        Y <- Y/matrix(sY, nrow(Y), p)
        w <- w * sY
		if (use.shortcut) maxit <- 0
		if (missing(MaxNWts)) MaxNWts <- (r + 1) * p
        if (length(offset) > 1L) {
            if (ncol(offset) != p)
                stop("ncol(offset) is wrong")
            mask <- c(rep(FALSE, r + 1L + p), rep(c(FALSE, rep(TRUE,
                r), rep(FALSE, p)), p - 1L))
            X <- cbind(X, offset)
            Wts <- as.vector(rbind(matrix(0, r + 1L, p), diag(p)))
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask,
                size = 0, skip = TRUE, softmax = TRUE, censored = FALSE,
                rang = 0, maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
        }
        else {
            mask <- c(rep(FALSE, r + 1L), rep(c(FALSE, rep(TRUE,
                r)), p - 1L))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0,
                skip = TRUE, softmax = TRUE, censored = FALSE,
                rang = 0, maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
        }
        if (use.shortcut) {
 			if (Xr==1) {
				if (missing(weights)) fit$fitted.values  <- matrix(apply(Y, 2, mean), nrow=nrow(Y),
				ncol=ncol(Y), byrow=TRUE, dimnames=dimnames(Y))
				else fit$fitted.values  <- matrix(apply(Y, 2, weighted.mean, w = c(w)), nrow=nrow(Y),
				ncol=ncol(Y), byrow=TRUE, dimnames=dimnames(Y))
			}
			else if (Xr==nrow(Y)) fit$fitted.values  <-  Y
 			if (missing(weights)) fit$value <- sum(-fit$fitted.values * log(fit$fitted.values), na.rm = TRUE)
 			else fit$value <- sum(-fit$fitted.values * log(fit$fitted.values) * c(w), na.rm = TRUE)
 		}
    }
    else {
        if (length(offset) <= 1L) {
            mask <- c(FALSE, rep(TRUE, r))
            fit <- nnet.default(X, Y, w, mask = mask, size = 0,
				skip = TRUE, entropy = TRUE, rang = 0, maxit = maxit,
				MaxNWts = MaxNWts, ..., PACKAGE="nnet")
        }
        else {
            mask <- c(FALSE, rep(TRUE, r), FALSE)
            Wts <- c(rep(0, r + 1L), 1)
            X <- cbind(X, offset)
            fit <- nnet.default(X, Y, w, Wts = Wts, mask = mask,
				size = 0, skip = TRUE, entropy = TRUE, rang = 0,
				maxit = maxit, MaxNWts = MaxNWts, ..., PACKAGE="nnet")
        }
    }
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$terms <- Terms
    fit$call <- call
    fit$weights <- w
    fit$deviance <- 2 * fit$value
	fit$entropy <-  fit$deviance/2/sum(c(w))
	fit$diversity <- exp(fit$entropy)
    fit$rank <- Xr
    edf <- (ncol(Y) - 1) * Xr
    if (length(dn <- colnames(Y)) > 0)
        fit$lab <- dn
    else fit$lab <- 1L:ncol(Y)
    fit$coefnames <- colnames(X)
    fit$vcoefnames <- fit$coefnames[1L:r]
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    fit$edf <- edf
    fit$AIC <- fit$deviance + 2 * edf
    if (model) fit$model <- m
    class(fit) <- c("mdm", "multinom", "nnet")
    fit
}

