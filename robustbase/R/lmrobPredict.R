#  File .../lmrobPredict.R
#  Part of the R package 'robustbase', http://www.R-project.org
#  Based on predict.lm (cf. src/library/stats/R/lm.R)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# Note that '# *rob' indicate adjustment for the robust case

predict.lmrob <-
    function(object, newdata, se.fit = FALSE, scale = NULL, df = NULL,  # *rob
	     interval = c("none", "confidence", "prediction"),
	     level = .95,  type = c("response", "terms"),
	     terms = NULL, na.action = na.pass, pred.var = res.var/weights,
             weights = 1, ...)
{
    tt <- terms(object)
    if(!inherits(object, "lmrob") && !inherits(object, "glmrob")) # *rob
	warning("calling predict.lm(<fake-lmrob-object>) ...") # *rob
    if(missing(newdata) || is.null(newdata)) {
	mm <- X <- model.matrix.lm(object)
	mmDone <- TRUE
	offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
                         xlev = object$xlevels)
        if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        offset <- rep.int(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset")))
            for(i in off.num)
                offset <- offset + eval(attr(tt, "variables")[[i+1]], newdata)
	if (!is.null(object$call$offset))
	    offset <- offset + eval(object$call$offset, newdata)
	mmDone <- FALSE
    }
    n <- length(object$residuals) # NROW(qr(object)$qr)
    p <- object$rank
    if(is.null(p)) { # *rob
        df <- Inf
        p <- sum(!is.na(coef(object)))
        piv <- seq_len(p)
    } else {
        p1 <- seq_len(p)
        piv <- if(p) qr(object)$pivot[p1]
    }
    if(p < ncol(X) && !(missing(newdata) || is.null(newdata)))
	warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    X.piv <- X[, piv, drop = FALSE]
    predictor <- drop(X.piv %*% beta[piv])
    if (!is.null(offset))
	predictor <- predictor + offset

    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata)) { # *rob: this and next if statement are combined
            warning("Predictions on current data refer to _future_ responses")
            if (missing(weights)) {
                w <- weights(object) # *rob
                if (!is.null(w)) {
                    weights <- w
                    warning("Assuming prediction variance inversely proportional to weights used for fitting")
                }
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) && missing(pred.var))
            warning("Assuming constant prediction variance even though model fit is weighted")
        if (inherits(weights, "formula")){
            if (length(weights) != 2L)
                stop("'weights' as formula should be one-sided")
            d <- if(missing(newdata) || is.null(newdata))
                model.frame(object)
            else
                newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }## "prediction" interval

    type <- match.arg(type)
    if(se.fit || interval != "none") {# *rob: whole 'then' statement is different
        df <- object$df.residual
	res.var <- if (is.null(scale)) object$s^2  else scale^2
	ip <- if(type != "terms")
	    diag(X.piv %*% object$cov %*% t(X.piv)) else rep.int(0, n)
    }

    if (type == "terms") { ## type == "terms" ------------

	if(!mmDone){
            mm <- model.matrix.lm(object) # *rob: call of model.matrix.lm
                                        # instead of model.matrix
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
	ll <- attr(tt, "term.labels")
	hasintercept <- attr(tt, "intercept") > 0L
	if (hasintercept) ll <- c("(Intercept)", ll)
	aaa <- factor(aa, labels = ll)
	asgn <- split(order(aa), aaa)
	if (hasintercept) {
	    asgn$"(Intercept)" <- NULL
	    if(!mmDone){
                mm <- model.matrix.lm(object) # *rob: call of model.matrix.lm
                                        # instead of model.matrix
                mmDone <- TRUE
            }
	    avx <- colMeans(mm)
	    termsconst <- sum(avx[piv] * beta[piv])
	}
	nterms <- length(asgn)
        if(nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))

            if (se.fit || interval != "none") {
                ip <- predictor # *rob: just this assignment is needed
            }
             if(hasintercept)
                X <- sweep(X, 2L, avx, check.margin=FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]      # Columns of X, ith term
                ii <- unpiv[iipiv]      # Corresponding rows of cov
                iipiv[ii == 0L] <- 0L
                predictor[, i] <-
                    if(any(iipiv > 0L)) X[, iipiv, drop = FALSE] %*% beta[iipiv]
                    else 0
                if (se.fit || interval != "none"){
                    ip[, i] <- if(any(iipiv > 0L)){# *rob: next steps modified
                        h.X <- X[, iipiv, drop = FALSE]
                        diag(h.X %*% object$cov[ii, ii] %*% t(h.X))
                    } else 0
                }
                }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit)
                    ip <- ip[, terms, drop = FALSE]
            }
        } else {                        # no terms
            predictor <- ip <- matrix(0, n, 0L)
        }
	attr(predictor, 'constant') <- if (hasintercept) termsconst else 0
    }

### Now construct elements of the list that will be returned

    if(interval != "none") {
	tfrac <- qt((1 - level)/2, df)
	hwid <- tfrac * switch(interval,
			       confidence = sqrt(ip),
			       prediction = sqrt(ip+pred.var)
			       )
	if(type != "terms") {
	    predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
	    colnames(predictor) <- c("fit", "lwr", "upr")
	} else {
            if (!is.null(terms)) hwid <- hwid[, terms, drop = FALSE]
	    lwr <- predictor + hwid
	    upr <- predictor - hwid
	}
    }
    if(se.fit || interval != "none") {
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms)) se <- se[, terms, drop = FALSE]
    }
    if(missing(newdata) && !is.null(na.act <- object$na.action)) {
	predictor <- napredict(na.act, predictor)
	if(se.fit) se <- napredict(na.act, se)
    }
    if(type == "terms" && interval != "none") {
	if(missing(newdata) && !is.null(na.act)) {
	    lwr <- napredict(na.act, lwr)
	    upr <- napredict(na.act, upr)
	}
	list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
	     df = df, residual.scale = sqrt(res.var))
    } else if (se.fit)
        list(fit = predictor, se.fit = se,
             df = df, residual.scale = sqrt(res.var))
    else predictor
}


