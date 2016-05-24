`get.models` <-
function(object, subset, cluster = NA, ...) {
    if (!inherits(object, "model.selection"))
		stop("'object' must be a 'model.selection' object")

	hasModelList <- is.list(attr(object, "modelList"))
	calls <- attr(object, "model.calls")
	if((hasNoCalls <- is.null(calls)) && !hasModelList)
		stop("'object' has no 'model.calls' attribute")

	if(!missing(subset)) {
		r <- subset_eval(substitute(subset), object, parent.frame())

		if(!isTRUE(r) && !is.na(r)) {
			if(is.character(r)) r <- match(r, dimnames(object)[[1L]])
		} else r <- TRUE
	} else {
		stop("'subset' is missing (use subset=TRUE to evaluate all models)")
	}

	newargs <- match.call()
	newargs[[1L]] <- NULL
	newargs[c('object', 'subset', 'cluster')] <- NULL
	naNames <- names(newargs)

	if(hasModelList) {
		DebugPrint(hasModelList)
		if(length(newargs) == 0L) {
			models <- attr(object, "modelList")[r]
			attr(models, "rank") <- attr(object, "rank")
			return(models)
		}
		DebugPrint("refitting...")
		if(hasNoCalls) calls <- lapply(attr(object, "modelList")[r], get_call)
	} else calls <- calls[r]

	if(length(newargs)) for(i in seq_along(calls)) calls[[i]][naNames] <- newargs

	doParallel <- inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck()
		# all this is to trick the R-check
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		models <- clusterApply(cluster, calls, "eval", envir = .GlobalEnv)
	} else {
	glo <- attr(object, "global")
	if(is.null(glo)) {
		models <- lapply(attr(object, "model.calls"), function(cl) {
			eval(cl, envir = environment(formula(cl)))
		})
	} else {
		env <- attr(tryCatch(terms(glo), error = function(...) terms(formula(glo))),
			".Environment")
		models <- lapply(calls, eval, envir = env)
	}
	}

	for(i in c("rank", "beta"))
		attr(models, i) <- attr(object, i)
	return(models)
}

`pget.models` <-
function(object, cluster = NA, subset, ...) {
	.Deprecated("get.models")
	get.models(object, subset, cluster, ...)
}
