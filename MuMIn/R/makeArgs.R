# 'makeArgs' internal generic to generate arguments from combination of
# term names (a character vector), additional list of arbitrary options is accepted.
# This is much a reverse action to getAllTerms
makeArgs <- function(obj, termNames, opt, ...) UseMethod("makeArgs", obj)

# opt == argsOptions

#argsOptions <- list(
#	response = attr(allTerms0, "response"),
#	intercept = nInts,        ### ONLY .default
#	interceptLabel = interceptLabel,
#	random = attr(allTerms0, "random"),
#	gmCall = gmCall,          ### ONLY .default
#	gmEnv = gmEnv,
#	allTerms = allTerms0,
#	gmCoefNames = gmCoefNames,
#	gmDataHead = if(!is.null(gmCall$data)) {
#		if(eval(call("is.data.frame", gmCall$data), gmEnv))
#			eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
#		} else NULL,
#	gmFormulaEnv = gmFormulaEnv
#	)


.getCoefNames <- 
function(formula, data, contrasts, envir = parent.frame()) {
	colnames(eval(call("model.matrix.default",
		object = formula, data = data, contrasts.arg = contrasts), envir = envir))
}

makeArgs.default <- 
function(obj, termNames, opt, ...) {
	reportProblems <- character(0L)
	termNames[termNames %in% opt$interceptLabel] <- "1"
	## XXX: what if length(opt$intercept) > 1 ???
	f <- reformulate(c(if(!opt$intercept) "0", termNames), response = opt$response)
	environment(f) <- opt$gmFormulaEnv
	ret <- list(formula = f)
	if(!is.null(opt$gmCall$start)) {
		coefNames <- fixCoefNames(.getCoefNames(f, opt$gmDataHead,
			opt$gmCall$contrasts, envir = opt$gmEnv))
		idx <- match(coefNames, opt$gmCoefNames)
		if(anyNA(idx)) reportProblems <-
			append(reportProblems, "cannot subset 'start' argument. Coefficients in model do not exist in 'global.model'")
		else ret$start <- substitute(start[idx], list(start = opt$gmCall$start,
			idx = idx))
	}
	#attr(ret, "formulaList") <- list(f)
	attr(ret, "problems") <- reportProblems
	ret
}

makeArgs.gls <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	names(ret)[1L] <- "model"
	ret
}

`makeArgs.asreml` <- 
makeArgs.MCMCglmm <-
makeArgs.lme <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	names(ret)[1L] <- "fixed"
	ret
}


`makeArgs.glmmadmb` <- 
`makeArgs.clmm` <- 		## Class 'clmm'  from package 'ordinal':
`makeArgs.merMod` <-    ## since lme4-0.99999911-0
`makeArgs.mer` <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	if(!is.null(opt$random)) ret[['formula']] <-
		update.formula(ret[['formula']], opt$random)
	ret
}


# used by makeArgs.unmarkedFit*
`.makeUnmarkedFitFunnyFormulas` <- 
function(termNames, opt, fnames) {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Int)", "(1)", termNames[i], fixed = TRUE)
	zexpr <- lapply(termNames, function(x) parse(text = x)[[1L]])
	zsplt <- split(sapply(zexpr, "[[", 2L), as.character(sapply(zexpr, "[[", 1L)))
	zarg <- lapply(zsplt, function(x) reformulate(as.character(x)))
	zarg <- zarg[fnames]
	names(zarg) <- fnames
	i <- sapply(zarg, is.null)
	zarg[i] <- rep(list(~ 1), sum(i))
	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	#print(zarg)
	zarg
}

# This works for all the child classes provided that the arguments containing
# formulas are named *formula, and they are at the beginnning of the call, and
# if there is no more than 6 of them, and if there is a single 'formula'
# argument it is doubly-one-sided.
# Guessing the argument names is only negligibly slower than when they are
# provided to the function, so the methods for child classes are commented out.

`makeArgs.unmarkedFit` <- 
function(obj, termNames, opt,
	fNames = sapply(obj@estimates@estimates, slot, "short.name"),
	formula.arg = "formula",
	...) {
	zarg <- .makeUnmarkedFitFunnyFormulas(termNames, opt, fNames)
	if(missing(formula.arg)) { # fallback
		nm <- names(opt$gmCall)[-1L]
		formula.arg <- nm[grep(".*formula$", nm[1L:7L])]
		#print(formula.arg)
	}
	if(length(formula.arg) == 1L && formula.arg == "formula") {
		i <- sapply(zarg, is.null)
		zarg[i] <- rep(list(~ 1), sum(i))
		ret <- list(formula = call("~", zarg[[2L]], zarg[[1L]][[2L]]))
		#attr(ret, "formulaList") <- zarg
		ret
	} else {
		names(zarg) <- formula.arg
		zarg
	}
}

`makeArgs.unmarkedFitDS` <-
function(obj, termNames, opt, ...)  {
	termNames <- sub("^p\\(sigma(.+)\\)", "p(\\1)", termNames, perl = TRUE)
	termNames[termNames == "p((Intercept))"] <- "p(1)"
	makeArgs.unmarkedFit(obj, termNames, opt, c("lam", "p"),
		"formula")
}

`makeArgs.coxph` <- 
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	ret$formula <- update.formula(ret$formula, . ~ . + 1)
	ret
}

`makeArgs.betareg` <- 
function(obj, termNames, opt, ...) {
	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Intercept)", "1", termNames[i], fixed = TRUE)
	j <- grepl("^\\(phi\\)_", termNames)
	zarg <- list(	
		beta = formula(terms.formula(reformulate(termNames[!j]), simplify = TRUE))
		)
	if(any(j))
		zarg$phi <- 
			formula(terms.formula(reformulate(substring(termNames[j], 7L)), simplify = TRUE))
	
	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	fexpl <- zarg$beta[[2L]]
	if(!is.null(zarg$phi)) fexpl <- call("|", fexpl, zarg$phi[[2L]])
		else zarg$phi <- NULL
	ret <- list(formula = call("~", opt$response, fexpl))
	#attr(ret, "formulaList") <- zarg
	ret
}

`makeArgs.hurdle` <- 
`makeArgs.zeroinfl` <-
function(obj, termNames, opt, ...) {

	intType <- substring(opt$interceptLabel, 0,
		regexpr("_", opt$interceptLabel, fixed = TRUE) - 1L)

	i <- termNames %in% opt$interceptLabel
	termNames[i] <- gsub("(Intercept)", "1", termNames[i], fixed = TRUE)
	pos <- regexpr("_", termNames, fixed = TRUE)

	fnames <- c("count", "zero")
	
	zarg <- split(substring(termNames, pos + 1L, 256L),
		substring(termNames, 1L, pos - 1L))
	for(j in fnames) zarg[[j]] <-
		if(is.null(zarg[[j]])) {
			if(j %in% intType) ~1 else ~0
		} else formula(terms.formula(reformulate(as.character(zarg[[j]]),
			intercept = j %in% intType),
			simplify = TRUE))

	zarg <- lapply(zarg, `environment<-`, opt$gmFormulaEnv)
	zarg <- zarg[fnames]
	fexpl <- zarg$count[[2L]]
	if(!is.null(zarg$zero)) fexpl <-
		call("|", fexpl, zarg$zero[[2L]]) else zarg$zero <- NULL
	ret <- list(formula = call("~", opt$response, fexpl))
	#attr(ret, "formulaList") <- zarg
	ret
}

#`makeArgs.unmarkedFitColExt` <- function(obj, termNames, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, opt, c("psi", "col", "ext", "p"),
#		c("psiformula", "gammaformula", "epsilonformula", "pformula"))
#
#
#`makeArgs.unmarkedFitGMM` <- function(obj, termNames, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, opt,
#		c("lambda", "phi", "p"),
#		c("lambdaformula", "phiformula", "pformula") )
#
#`makeArgs.unmarkedFitPCO` <- function(obj, termNames, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, opt,
#		c("lam", "gamConst", "omega", "p"),
#		c("lambdaformula", "gammaformula", "omegaformula", "pformula"))
#
#`makeArgs.unmarkedFitOccu` <- function(obj, termNames, opt, ...)
#	makeArgs.unmarkedFit(obj, termNames, opt, c("psi", "p"),
#		"formula")


`makeArgs.coxme` <-
`makeArgs.lmekin` <-
function(obj, termNames, opt, ...) {
	ret <- makeArgs.default(obj, termNames, opt)
	ret$formula <- update.formula(update.formula(ret$formula, . ~ . + 1),
		opt$random)
	ret
}

`makeArgs.mark` <- 
function(obj, termNames, opt, ...) {
	interceptLabel <- "(Intercept)"
	termNames <- sub(interceptLabel, "1", termNames, fixed = TRUE)
	rxres <- regexpr("^([a-zA-Z]+)\\((.*)\\)$", termNames, perl = TRUE)
	cs <- attr(rxres, "capture.start")
	cl <- attr(rxres, "capture.length")
	parname <- substring(termNames, cs[, 1L], cs[, 1L] + cl[,1L] - 1L)
	parval <- substring(termNames, cs[, 2L], cs[, 2L] + cl[,2L] - 1L)
	
	#formulaList <- lapply(split(parval, parname), function(x) {
	#	int <- x == "1"
	#	x <- x[!int]
	#	res <- if(!length(x))
	#			if(int) ~ 1 else ~ 0 else 
	#		reformulate(x, intercept = any(int))
	#	environment(res) <- opt$gmFormulaEnv
	#	res
	#})
	
	mpar <- if(is.null(obj$model.parameters))
		eval(opt$gmCall$model$parameters) else
		obj$model.parameters
	#for(i in names(mpar)) mpar[[i]]$formula <- formulaList[[i]]
	#ret <- list(model.parameters = mpar)
	
	if(opt$gmCall[[1L]] == "run.mark.model") {
		arg.model <- opt$gmCall$model
		arg.model$parameters <- mpar
		ret <- list(model = arg.model)
	} else {
		ret <- list(model.parameters = mpar)
	}
	
	#attr(ret, "formulaList") <- formulaList
	ret
}


`makeArgs.aodml` <-
function(obj, termNames, opt, ...) {
	if(sys.nframe() > 2L && (parent.call <- sys.call(-2L))[[1L]] == "dredge" &&
	   !is.null(get_call(obj)$fixpar))
		stop(simpleError("'aodml' models with constant parameters cannot be handled by 'dredge'",
						 call = parent.call))
	makeArgs.default(obj, termNames, opt, ...)
}



