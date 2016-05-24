`model.avg` <-
function (object, ..., revised.var = TRUE) {
	if (isTRUE("method" %in% names(match.call())))
		stop("argument 'method' is defunct")
	UseMethod("model.avg")
}

.coefarr.avg <-
function(cfarr, weight, revised.var, full, alpha) {	
	weight <- weight / sum(weight)
	nCoef <- dim(cfarr)[3L]
	if(full) {
		nas <- is.na(cfarr[, 1L, ]) & is.na(cfarr[, 2L, ])
		cfarr[, 1L, ][nas] <- cfarr[, 2L, ][nas] <- 0
		#cfarr[, 1L:2L, ][is.na(cfarr[, 1L:2L, ])] <- 0
		if(!all(is.na(cfarr[, 3L, ])))
			cfarr[ ,3L, ][is.na(cfarr[ , 3L, ])] <- Inf
	}
	
	avgcoef <- array(dim = c(nCoef, 5L),
		dimnames = list(dimnames(cfarr)[[3L]], c("Estimate",
		"Std. Error", "Adjusted SE", "Lower CI", "Upper CI")))
	for(i in seq_len(nCoef))
		avgcoef[i, ] <- par.avg(cfarr[, 1L, i], cfarr[, 2L, i], weight,
			df = cfarr[, 3L, i], alpha = alpha, revised.var = revised.var)
		
	avgcoef[is.nan(avgcoef)] <- NA
	return(avgcoef)
}

`model.avg.model.selection` <-
function(object, subset, fit = FALSE, ..., revised.var = TRUE) {

	if(!missing(subset)) {
		cl <- match.call()
		cl[[1L]] <- as.name("subset")
		names(cl)[2L] <- "x"
		object <- eval.parent(cl[1L:3L])
	}
	
	# TODO: unify refitting conditions in model.avg and model.sel
	
	if(fit || !missing(...)) {
		cl <- match.call()
		cl$fit <- NULL
		arg1 <- names(cl)[-(1L:2L)] %in% names(formals("model.avg.default"))
		cl1 <- cl[c(TRUE, TRUE, !arg1)]
		cl1[[1L]] <- as.name("get.models")
		if(is.null(cl1[["subset"]])) cl1[["subset"]] <- NA
		# TODO: subset = TRUE

		cl2 <- cl[c(TRUE, TRUE, arg1)]
		cl2[[2L]] <- cl1
		cl2[[1L]] <- as.name("model.avg")
		#message("Re-fitting model objects...")
		return(eval(cl2, parent.frame()))
	}

	if(nrow(object) <= 1L) stop("'object' consists of only one model")
	
	ct <- attr(object, "coefTables")

	cfarr <- coefArray(ct)
	weight <- Weights(object)

	cfmat <- cfarr[, 1L, ]
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), colnames(cfmat)))
	coefMat[1L, ] <- drop(weight %*% cfmat)
	coefMat[2L, ] <- coefMat[1L, ] / colSums(array(weight *
		as.numeric(!is.na(cfarr[, 1L, ])), dim = dim(cfmat)))
	coefMat[is.nan(coefMat)] <- NA_real_

	#allterms1 <- lapply(attr(object, "calls"), function(x)
		#getAllTerms(as.formula(x[[switch(as.character(x[[1L]]),
			#lme=, lme.formula= "fixed", gls= "model", "formula")]])))
	all.terms <- attr(object, "terms")
	all.vterms <- all.terms[!(all.terms %in% attr(all.terms, "interceptLabel")
		| apply(is.na(object[, all.terms]), 2L, all))]
	#allterms1 <- apply(!is.na(object[, all.vterms, drop = FALSE]), 1L, function(x) all.vterms[x])
	allterms1 <- applyrns(!is.na(object[, all.vterms, drop = FALSE]), function(x) all.vterms[x])
	allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.vterms)

	mstab <- itemByType(object, c("df", "loglik", "ic", "delta", "weight"))
	rownames(mstab) <- allmodelnames
	

	ret <- list(
		msTable = structure(as.data.frame(mstab),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = cfarr,
		importance = importance(object),
		x = NULL,
		residuals = NULL,
		formula = if(!is.null(attr(object, "global")))
			formula(attr(object, "global")) else NULL,
		call = {
			cl <- match.call()
			cl[[1L]] <- as.name(.Generic)
			cl
		}
	)

	attr(ret, "rank") <- attr(object, "rank")
	attr(ret, "modelList") <- attr(object, "modelList")
	if(is.null(attr(ret, "modelList")))
		attr(ret, "model.calls") <- attr(object, "model.calls")
	attr(ret, "beta") <- attr(object, "beta")
	attr(ret, "nobs") <- attr(object, "nobs")
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

`model.avg.default` <-
function(object, ..., beta = c("none", "sd", "partial.sd"),
		 rank = NULL, rank.args = NULL, revised.var = TRUE,
		 dispersion = NULL, ct.args = NULL) {

	if (is.object(object)) {
		models <- list(object, ...)
        rank <- .getRank(rank, rank.args = rank.args, object = object) 
	} else {
		if(length(object) == 0L) stop("'object' is an empty list")
		models <- object
		object <- object[[1L]]
        if (!is.null(rank) || is.null(rank <- attr(models, "rank"))) {
            rank <- .getRank(rank, rank.args = rank.args, object = object)
      	}
	}
	
	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)

	nModels <- length(models)
	if(nModels == 1L) stop("only one model supplied. Nothing to do")
	.checkModels(models)
	
	testSmoothKConsistency(models) # for gam, if any

	ICname <- asChar(attr(rank, "call")[[1L]])

	allterms1 <- lapply(models, getAllTerms)
	all.terms <- unique(unlist(allterms1, use.names = FALSE))

	# sort by level (main effects first)
	all.terms <- all.terms[order(vapply(gregexpr(":", all.terms),
		function(x) if(x[1L] == -1L) 0L else length(x), 1L), all.terms)]

	# allmodelnames <- modelNames(models, asNumeric = FALSE,
		# withRandomTerms = FALSE, withFamily = FALSE)
	allmodelnames <- .modelNames(allTerms = allterms1, uqTerms = all.terms)

	#if(is.null(names(models))) names(models) <- allmodelnames
	
	coefTableCall <- as.call(c(alist(coefTable, models[[i]],
		dispersion = dispersion[i]), ct.args))
	if(is.null(dispersion)) coefTableCall$dispersion <- NULL
	
	if(betaMode == 2L) {
		coefTableCall[[1L]] <- as.name("std.coef")
		coefTableCall[['partial.sd']] <- TRUE
	}
	
	DebugPrint(coefTableCall)

	# check if models are unique:
	if(!is.null(dispersion)) dispersion <- rep(dispersion, length.out = nModels)
	coefTables <- vector(nModels, mode = "list")
	for(i in seq_len(nModels))
		coefTables[[i]] <-  eval(coefTableCall)
	
	mcoeffs <- lapply(coefTables, "[", , 1L)
	dup <- unique(sapply(mcoeffs, function(i) which(sapply(mcoeffs, identical, i))))
	dup <- dup[sapply(dup, length) > 1L]
	if (length(dup) > 0L) stop("models are not unique. Duplicates: ",
		prettyEnumStr(sapply(dup, paste0, collapse = " = "),
			quote = "'"))

	LL <- .getLik(object)
	logLik <- LL$logLik
	lLName <- LL$name

	ic <- vapply(models, rank, 0)
	logLiks <- lapply(models, logLik)
	delta <- ic - min(ic)
	weight <- exp(-delta / 2) / sum(exp(-delta / 2))
	model.order <- order(weight, decreasing = TRUE)

	# ----!!! From now on, everything MUST BE ORDERED by 'weight' !!!-----------
	mstab <- cbind(df = vapply(logLiks, attr, 0, "df"),
		logLik = as.numeric(logLiks), IC = ic, delta = delta, weight = weight,
		deparse.level = 0L)
	if(!is.null(dispersion)) mstab <- cbind(mstab, Dispersion = dispersion)
	rownames(mstab) <- allmodelnames
	mstab <- mstab[model.order, ]
	weight <- mstab[, "weight"] # has been sorted in table
	models <- models[model.order]
	coefTables <- coefTables[model.order]

	if (betaMode == 1L) {
		response.sd <- sd(model.response(model.frame(object)))
		for(i in seq_along(coefTables))
			coefTables[[i]][, 1L:2L] <-
				coefTables[[i]][, 1L:2L] *
				apply(model.matrix(models[[i]]), 2L, sd) / response.sd
	}

	cfarr <- coefArray(coefTables)
	cfmat <- cfarr[, 1L, ]
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), colnames(cfmat)))
	coefMat[1L, ] <- drop(weight %*% cfmat)
	coefMat[2L, ] <- coefMat[1L, ] / colSums(array(weight *
		as.numeric(!is.na(cfarr[, 1L, ])), dim = dim(cfmat)))
	coefMat[is.nan(coefMat)] <- NA_real_
		

    names(all.terms) <- seq_along(all.terms)
	colnames(mstab)[3L] <- ICname

	# Benchmark: 3.7x faster
	#system.time(for(i in 1:10000) t(array(unlist(p), dim=c(length(all.terms),length(models)))))
	#system.time(for(i in 1:10000) do.call("rbind", p))

	vpresent <- do.call("rbind", lapply(models, function(x)
		all.terms %in% getAllTerms(x)))

	importance <- apply(weight * vpresent, 2L, sum)
	names(importance) <- all.terms
	o <- order(importance, decreasing = TRUE)
	importance <- importance[o]
	attr(importance, "n.models") <- structure(colSums(vpresent)[o], names = all.terms)
	class(importance) <- c("importance", "numeric") 

	mmxs <- tryCatch(cbindDataFrameList(lapply(models, model.matrix)),
					 error = return_null, warning = return_null)

	# Far less efficient:
	#mmxs <- lapply(models, model.matrix)
	#mx <- mmxs[[1]];
	#for (i in mmxs[-1])
	#	mx <- cbind(mx, i[,!(colnames(i) %in% colnames(mx)), drop=FALSE])

	# residuals averaged (with brute force)
	#rsd <- tryCatch(apply(vapply(models, residuals, residuals(object)), 1L,
		#weighted.mean, w = weight), error = return_null)
	#rsd <- NULL
	## XXX: how to calc residuals ?
	
	modelClasses <- lapply(models, class)
	frm <-
	if(all(vapply(modelClasses[-1L], identical, FALSE, modelClasses[[1L]]))) {
		trm <- tryCatch(terms(models[[1L]]),
				error = function(e) terms(formula(models[[1L]])))
		response <- attr(trm, "response")
		m1 <- models[[1L]]
		makeArgs(m1, all.terms, opt = list(
			response = if(response > 0L) attr(trm, "variables")[[response + 1L]] else NULL,
			gmFormulaEnv = environment(formula(m1)),
			intercept = ! identical(unique(unlist(lapply(allterms1, attr, "intercept"))), 0),
			interceptLabel = unique(unlist(lapply(allterms1, attr, "interceptLabel"))),
			#	random = attr(allTerms0, "random"),
			gmCall = get_call(m1),
			gmEnv = parent.frame(),
			allTerms = all.terms,
			random = . ~ .
			))[[1L]]
	} else NA
	

	ret <- list(
		msTable = structure(as.data.frame(mstab),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = cfarr,
		importance = importance,
		x = mmxs,
		residuals = NULL, # no residuals
		formula = frm,
		call = {
			cl <- match.call()
			cl[[1L]] <- as.name(.Generic)
			cl
		}
	)

	attr(ret, "rank") <- rank
	attr(ret, "modelList") <- models
	attr(ret, "beta") <- beta
	attr(ret, "nobs") <- nobs(object)
	attr(ret, "revised.var") <- revised.var
	class(ret) <- "averaging"
	return(ret)
}

.checkFull <- function(object, full, warn = TRUE) {
	if(isTRUE(attr(object, "arm")) && !full) {
		if(warn) cry(-1L, "'subset' averaged coefficients are not available with ARM algorithm",
					 warn = TRUE)
		return(TRUE)
	} else return(full)
}

`coef.averaging` <-
function(object, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 
	object$coefficients[if(full) 1L else 2L, ]
}

`predict.averaging` <-
function(object, newdata = NULL, se.fit = FALSE, interval = NULL,
	type = NA, backtransform = FALSE,
	full = TRUE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 


	if (!missing(interval)) .NotYetUsed("interval", error = FALSE)
	
	if(backtransform && !is.na(type) && type == "response")
		warning("back-transforming predictions that are already on response scale")

	models <- attr(object, "modelList")
	if(is.null(models)) stop("can predict only from 'averaging' object containing model list")

	# If all models inherit from lm:
	if ((missing(se.fit) || !se.fit)
		&& (is.na(type) || type == "link")
		&& !backtransform
		&& all(linherits(models, c(gam = FALSE, lm = TRUE)))
		&& !anyNA(object$coefficients[1L, ])
		) {
		
		coeff <- coef(object, full = full)

		X <- model.matrix(object)
		if (missing(newdata) || is.null(newdata)) {
			Xnew <- X
		} else {
			tt <- delete.response(terms(formula(object)))
			xlev <- unlist(unname(lapply(models, "[[", "xlevels")),
						   recursive = FALSE, use.names = TRUE)
			Xnew <- model.matrix(tt, data = newdata, xlev = xlev)
		}
		
		colnames(Xnew) <- fixCoefNames(colnames(Xnew))
		Xnew <- Xnew[, match(names(coeff), colnames(Xnew), nomatch = 0L)]
		ret <- (Xnew %*% coeff)[, 1L]
		#if (se.fit) {
		#	scale <- 1
		#	covmx <- solve(t(X) %*% X)
		#	se <- sqrt(diag(Xnew %*% covmx %*% t(Xnew))) * sqrt(scale) ## TODO: use matmult
		#	return(list(fit = y, se.fit = se))
		#}
	} else {
		# otherwise, use brute force:

		if(full == FALSE) warning("argument 'full' ignored")

		cl <- as.list(match.call())
		cl$backtransform <- cl$full <- NULL
		cl[[1L]] <- as.name("predict")
		cl <- as.call(cl)
		pred <- lapply(models, function(x, pfr) {
			cl[[2L]] <- x
			y <- tryCatch({
				y <- eval(cl, pfr)
				if(is.numeric(y)) y else structure(as.list(y[c(1L, 2L)]),
					names = c("fit", "se.fit"))
				}, error = function(e) e)
		}, pfr = parent.frame())
		
		err <- sapply(pred, inherits, "condition")
		if (any(err)) {
			lapply(pred[err], warning)
			stop(sprintf(ngettext(sum(err), "'predict' for model %s caused error",
				"'predict' for models %s caused errors"),
				prettyEnumStr(names(models[err]), quote = "'")
				))
		}

		.untransform <- function(fit, se.fit = NULL, models) {
			links <- tryCatch(vapply(models, function(m) family(m)[["link"]], ""),
							error = function(e) NULL)
			
			if (!is.null(links)) {
				if(any(links[1L] != links[-1L]))
					cry(sys.call(-1L), "cannot inverse-transform averaged prediction of models using different link functions")
				fam1 <- family(models[[1L]])
				if(is.null(se.fit)) return(fam1$linkinv(fit))
				else return(list(fit = fam1$linkinv(fit),
					se.fit = se.fit * abs(fam1$mu.eta(fit))
					))
			}
			return(NULL)
		}

		if (all(sapply(pred, is.list))) {
		#if(all(sapply(pred, function(x) c("fit", "se.fit") %in% names(x)))) {
			fit <- do.call("cbind", lapply(pred, "[[", "fit"))
			se.fit <- do.call("cbind", lapply(pred, "[[", "se.fit"))
			revised.var <- attr(object, "revised.var")

			apred <- unname(vapply(seq(nrow(fit)), function(i)
				par.avg(fit[i, ], se.fit[i, ], weight = Weights(object),
					df = NA_integer_, revised.var = revised.var),
					FUN.VALUE = double(5L)))

			# TODO: ase!
			#no.ase <- all(is.na(object$coefTable[,3]))
			# if(no.ase) 2 else 3
			ret <-  if (backtransform)
					.untransform(apred[1L, ], apred[2L, ], models = models) else
						list(fit = apred[1L, ], se.fit = apred[2L, ])
		} else {
			#tryCatch({
			i <- !vapply(pred, is.numeric, FALSE)
			if(any(i)) pred[i] <- lapply(pred[i], "[[", 1L)
			ret <- apply(do.call("cbind", pred), 1L, weighted.mean,
				w = Weights(object))
			if (backtransform)
				ret <- .untransform(ret, models = models)
		}

	}
	return(ret)
}

`fitted.averaging` <-
function (object, ...) .NotYetImplemented()
# predict.averaging(object, backtransform = TRUE, type = "link")

model.frame.averaging <-
function (formula, ...) {
	mergeMF(getModelList(formula))
}

model.matrix.averaging <-
function (object, ...) {
	if(j <- match("x", names(object), nomatch = 0L)) return(object[[j]])
	
	mf <- model.frame(object)
	do.call("model.matrix", list(object = terms(mf), data = mf,
		contrasts.arg = get.contrasts(object)))
}

`summary.averaging` <-
function (object, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)

	.makecoefmat <- function(cf) {
		no.ase <- all(is.na(cf[, 3L]))
		z <- abs(cf[, 1L] / cf[, if(no.ase) 2L else 3L])
		pval <- 2 * pnorm(z, lower.tail = FALSE)
		cbind(cf[, if(no.ase) 1L:2L else 1L:3L],
			`z value` = z, `Pr(>|z|)` = zapsmall(pval))
	}
	
	is.arm <- ncol(object$msTable) == 6L && (colnames(object$msTable)[6L] == "ARM weight")

	weight <- object$msTable[, if(is.arm) 6L else 5L]
	
	object$coefmat.full <- .makecoefmat(.coefarr.avg(object$coefArray, weight,
			attr(object, "revised.var"), TRUE, 0.05))

	if(!is.arm) object$coefmat.subset <-
		.makecoefmat(.coefarr.avg(object$coefArray, weight,
			attr(object, "revised.var"), FALSE, 0.05))
	
	object$coef.nmod <- colSums(!is.na(object$coefArray[, 1L,]))

	structure(object, ARM = is.arm, class = c("summary.averaging", "averaging"))
}

`confint.averaging` <-
function (object, parm, level = 0.95, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 


    a2 <- 1 - level
    a <- a2 / 2
    cf <- object$coefArray[, 1L, ]
    pnames <- colnames(cf)
    if (missing(parm)) parm <- pnames
		else if (is.numeric(parm)) parm <- pnames[parm]
	missing.par <- is.na(cf)
    se <- object$coefArray[, 2L, ]
    dfs <- object$coefArray[, 3L, ]
	if(full) {
		se[missing.par] <- cf[missing.par] <- 0
		if(!all(is.na(dfs))) dfs[missing.par] <- Inf
	}
    wts <- Weights(object) ## XXX: !
    ci <- t(sapply(parm, function(i)
		par.avg(cf[,i], se[,i], wts, dfs[, i], alpha = a2)))[, 4L:5L]
    
	
	ci[is.na(object$coefficients[1L, ]), ] <- NA_real_
    colnames(ci) <- getFrom("stats", "format.perc")(c(a, 1L - a), 3L)
    return(ci)
}

`print.summary.averaging` <-
function (x, digits = max(3L, getOption("digits") - 3L),
    signif.stars = getOption("show.signif.stars"), ...) {

    cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
		
	comcallstr <- 
	if(!is.null(attr(x, "model.calls"))) {
		commonCallStr(calls = attr(x, "model.calls"))
	} else if(!is.null(attr(x, "modelList"))) {
		commonCallStr(attr(x, "modelList"))
	} else NA
	if(!is.na(comcallstr)) {
		cat("Component model call: \n")
		cat(strwrap(comcallstr), sep = " \n     ")
	}		
		
    cat("\nComponent models: \n")
	print(round(as.matrix(x$msTable), 2L), na.print = "")

	if(!is.null(attr(x$msTable, "term.codes"))) {
		cat("\nTerm codes: \n")
		print.default(attr(x$msTable, "term.codes"), quote = FALSE)
	}
	cat("\nModel-averaged coefficients: ")
	if (nnotdef <- sum(is.na(x$coefmat.full[, 1L])))
		 cat("\n(", nnotdef, " not defined because of singularities in all ",
			"component models)", sep = "")

	hasPval <- TRUE
	coefTitles <- if(attr(x, "ARM"))
		c(coefmat.full = "(ARM average)") else
		c(coefmat.full = "(full average)",
		  coefmat.subset = "(conditional average)")
		
	n <- length(coefTitles)	
	for(i in seq.int(n)) {
		iname <- names(coefTitles[i])
		if(is.null(x[[i]])) next
		cat(" \n", coefTitles[i], " \n", sep = "")
		printCoefmat(x[[iname]], P.values = hasPval, has.Pvalue = hasPval,
			digits = digits, signif.stars = signif.stars,
			signif.legend = i == n)
	}
	
	#if (no.ase) cat("Confidence intervals are unadjusted \n")
	#printCoefmat(matrix(x$coef  .shrinkage, nrow = 1L,
		#dimnames = list("", x$term.names)), P.values = FALSE,
		#has.Pvalue = FALSE, cs.ind = seq_along(x$term.names), tst.ind = NULL)

	cat("\nRelative variable importance: \n")
	print(round(x$importance, 2L))
}

`print.averaging` <-
function(x, ...) {
	## XXX: backward compatibility:
	x <- upgrade_averaging_object(x)
    cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Component models:", "\n")
	comp.names <- rownames(x$msTable)
	comp.names[comp.names == ""] <- "null"
	cat(format(sQuote(comp.names), justify = "l"), fill = TRUE)
	cat("\nCoefficients:", "\n")
	print(x$coefficients[!is.na(x$coefficients[,1L]), , drop = FALSE])
    x
}

`vcov.averaging` <- function (object, full = FALSE, ...) {
	## XXX: backward compatibility:
	object <- upgrade_averaging_object(object)
	full <- .checkFull(object, full) 

	full <- as.logical(full)[1L]

	models <- attr(object, "modelList")
	if(is.null(models)) stop("cannot calculate covariance matrix from ",
							 "'averaging' object without component models")

	vcovs <- lapply(lapply(models, vcov), as.matrix)
	names.all <- dimnames(object$coefArray)[[3L]]
	nvars <- length(names.all)
	nvarseq <- seq(nvars)
	wts <- Weights(object)
	wts <- wts / sum(wts) # normalize just in case

	vcov0 <- matrix(if(full) 0 else NA_real_, nrow = nvars,
		ncol = nvars, dimnames = list(names.all, names.all))

	vcovs2 <- lapply(vcovs, function(v) {
		i <- match(fixCoefNames(dimnames(v)[[1L]]), names.all)
		vcov0[i, i] <- v
		return(vcov0)
	})
	b1 <- object$coefArray[, 1L, ]
	if(full) b1[is.na(b1)] <- 0
	avgb <- object$coefficients[2L - full, ]
	
	res <- sapply(nvarseq, function(c1) sapply(nvarseq, function(c2) {
		 weighted.mean(sapply(vcovs2, "[", c1, c2) + (b1[, c1] - avgb[c1]) *
			(b1[, c2] - avgb[c2]), wts, na.rm = TRUE)
	}))
	dimnames(res) <- list(names.all, names.all)
	return(res)
}

`logLik.averaging` <- function (object, ...) {
	models <- attr(object, "modelList")
	if(is.null(models)) {
		nobs <- attr(object, "nobs")
		apply(object$msTable, 1L, function(x) structure(list(x[2L]),
			df = x[1L], nobs = nobs, class = "logLik"))
	} else {
		structure(lapply(attr(object, "modelList"), logLik),
			names = rownames(object$msTable))
	}
}

`coefTable.averaging` <-
function (model, full = FALSE, adjust.se = TRUE, ...) {
	full <- .checkFull(model, full) 
	
	
    no.ase <- any(is.na(model$coefArray[,3L,]) & !is.na(model$coefArray[,1L,]))
	if(!missing(adjust.se) && adjust.se && no.ase) 
        warning("adjusted std. error not available for this type of model")
		
	weight <- model$msTable[, ncol(model$msTable)]

    cols <- c(1L, if(!adjust.se || no.ase) 2L else 3L)
	ct <- .coefarr.avg(model$coefArray, weight, TRUE, full, .05)[, cols] 
	.makeCoefTable(ct[,1L], ct[,2L], NA, rownames(ct))
}


## XXX: backward compatibility (< 0.15.0):
upgrade_averaging_object <-
function(x) {
	if(is.matrix(x$coefficients)) return(x)
	if(all(c("coefTable", "coef.shrinkage") %in% names(x))) {
		x$coefficients <- rbind(full = x$coef.shrinkage, subset = x$coefTable[, 1L])
		x$coefTable <- NULL
	} else stop("'averaging' object is corrupt")
	x
}
