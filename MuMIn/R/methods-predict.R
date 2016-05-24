## Predict methods for objects for which they are not available in their
## original packages, or replacements.

# Add 'se.fit' argument for predict:
# https://stat.ethz.ch/pipermail/r-help/2004-April/050144.html
# http://web.archiveorange.com/archive/v/rOz2zbtjRgntPMuIDoIl

# based on the original 'predict.gls' in package 'nlme'
`predict.gls` <-
function (object, newdata, se.fit = FALSE, na.action = na.fail, ...) {
    if (missing(newdata) && missing(se.fit)) return(fitted(object))
	
    form <- formula(object)[-2L]
	if(length(form[[2L]]) == 3L && form[[2L]][[1L]] == "|" )
		form[[2L]] <- form[[2L]][[2L]] 
	
	dataMod <- model.frame(object, data = newdata, na.action = na.action,
										drop.unused.levels = TRUE)
	contr <- object$contrasts
	for(i in names(contr)) {
		levs <- levels(dataMod[, i])
		if (any(wch <- is.na(match(levs, rownames(contr[[i]]))))) {
                stop(sprintf(ngettext(sum(wch), "level %s not allowed for %s", 
                  "levels %s not allowed for %s"), paste(levs[wch], 
                  collapse = ",")), domain = NA)
            }
		attr(dataMod[[i]], "contrasts") <- contr[[i]][levs, , drop = FALSE]
	}
	X <- model.matrix(terms(form), data = dataMod)
	
    cf <- coef(object)
    val <- c(X[, names(cf), drop = FALSE] %*% cf)
	if(se.fit) {
		# se <- sqrt(diag(X %*% vcov(object) %*% t(X)))
		# se <- sqrt(rowSums((X %*% vcov(object)) * X))
		se <- sqrt(matmultdiag(X %*% vcov(object), ty = X))
		val <- list(fit = val, se.fit = unname(se))
	}
	attr(val, "label") <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
        attr(val, "label") <- paste(attr(val, "label"), aux)
    }
	val
}

`predict.lme` <-
function (object, newdata, level, asList = FALSE,
	na.action = na.fail, se.fit = FALSE, ...) {
	cl <- match.call()
	cl$se.fit <- NULL
	cl[[1L]] <- call("get", "predict.lme", asNamespace("nlme"))	
	res <- eval.parent(cl)
	
	if(se.fit && (missing(level) || any(level > 0)))
		warning("cannot calculate standard errors for level > 0")
	if(se.fit && !missing(level) && length(level) == 1L && all(level == 0)) {
		if (missing(newdata) || is.null(newdata)) {
			X <- model.matrix(object, data = object$data)
		} else {
			tt <- delete.response(terms(formula(object)))
			xlev <- .getXlevels(tt, model.frame(object, data = object$data))
			X <- model.matrix(tt, data = newdata, contrasts.arg =
							  object$contrasts, xlev = xlev)
		}
		se <- sqrt(matmultdiag(X %*% vcov(object), ty = X))
		# se <- sqrt(rowSums((X %*% vcov(object)) * X))
		# se <- sqrt(diag(X %*% vcov(object) %*% t(X))) ## TODO: use matmult
		names(se) <- names(res)
		list(fit = c(res), se.fit = se)
	} else res
}

`predict.merMod` <-
function (object, newdata, type = c("link", "response"), se.fit = FALSE,
		re.form = NULL, ...) {
	
	if(!se.fit) {
		cl <- sys.call()
		cl$se.fit <- NULL
		return(do.call(getFrom("lme4", "predict.merMod"), as.list(cl[-1L]), envir = parent.frame()))
	}
	
	level0 <- (!is.null(re.form) && !inherits(re.form, "formula") && is.na(re.form)) || 
        (inherits(re.form, "formula") && length(re.form) == 2L && identical(re.form[[2L]], 
            0))
	if(!level0) stop("cannot calculate predictions with both standard errors and random effects")

	.predict_glm(object, newdata, type, se.fit,
			trms = delete.response(terms(formula(object, fixed.only = TRUE))),
			coeff = lme4::fixef(object),
			offset = lme4::getME(object, "offset"),
			...)
}


.predict_glm <-
function (object, newdata, type = c("link", "response"), se.fit = FALSE,
	trms, coeff, offset, ...) {
	
    type <- match.arg(type)
    if (!missing(newdata) && !is.null(newdata)) {
        xlev <- .getXlevels(trms, model.frame(trms, data = newdata))
        X <- model.matrix(trms, data = newdata, contrasts.arg = attr(model.matrix(object), 
            "contrasts"), xlev = xlev)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(trms, "offset"))) 
            for (i in off.num) offset <- offset + eval(attr(trms, 
                "variables")[[i + 1L]], newdata)
		
		cl <- get_call(object)
        if (!is.null(cl$offset)) 
            offset <- offset + eval(cl$offset, newdata)
    } else {
        X <- model.matrix(object)
        if (!length(offset))  offset <- NULL
    }
    y <- (X %*% coeff)[, 1L]
    if (!is.null(offset)) 
        y <- y + offset
    fam <- family(object)
    if (se.fit) {
        # covmat <- as.matrix(vcov(object))
		se <- sqrt(matmultdiag(X %*% as.matrix(vcov(object)), ty = X))
		# se <- sqrt(rowSums((X %*% covmat) * X))
        # se <- sqrt(diag(X %*% covmat %*% t(X)))
        if (type == "response" && inherits(fam, "family")) 
            list(fit = fam$linkinv(y), se.fit = se * abs(fam$mu.eta(y)))
        else list(fit = y, se.fit = se)
    } else {
        if (type == "response" && inherits(fam, "family")) 
            fam$linkinv(y)
        else y
    }
}

`predict.gamm` <- 
function (object, ...) mgcv::predict.gam(object[['gam']], ...)

