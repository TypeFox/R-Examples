## TODO: chunk size for evaluate = FALSE

`pdredge` <-
function(global.model, cluster = NA,
	beta = c("none", "sd", "partial.sd"),
	evaluate = TRUE,
	rank = "AICc", fixed = NULL, m.lim = NULL, m.min, m.max, subset,
	trace = FALSE, varying, extra, ct.args = NULL, check = FALSE, ...) {

#FIXME: m.max cannot be 0 - e.g. for intercept only model

	trace <- min(as.integer(trace), 2L)

	strbeta <- betaMode <- NULL
	eval(.expr_beta_arg)


###PAR
	qlen <- 25L
	# Imports: clusterCall, clusterApply
	doParallel <- evaluate && inherits(cluster, "cluster")
	if(doParallel) {
		.parallelPkgCheck() # XXX: workaround to avoid importing from 'parallel'
		clusterCall <- get("clusterCall")
		clusterApply <- get("clusterApply")
		clusterCall(cluster, "require", .packageName, character.only = TRUE)
		.getRow <- function(X) clusterApply(cluster, X, fun = ".pdredge_process_model")
	} else {
		.getRow <- function(X) lapply(X, pdredge_process_model, envir = props)
		clusterCall <- function(...) NULL
		message("Not using cluster.")
	}
###PAR

	gmEnv <- parent.frame()
	gmNobs <- nobs(global.model)

	gmCall <- get_call(global.model)
	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			stop("need a 'global.model' with a call component. Consider using ",
				if(inherits(global.model, c("gamm", "gamm4")))
					"'uGamm'" else "'updateable'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen
		if(!is.function(eval.parent(gmCall[[1L]])))
			cry(, "could not find function '%s'", asChar(gmCall[[1L]]))
	} else {
		# if 'update' method does not expand dots, we have a problem with
		# expressions like ..1, ..2 in the call. So try to replace them with
		# respective arguments in the original call
		isDotted <- grep("^\\.\\.", sapply(as.list(gmCall), asChar))
		if(length(isDotted) != 0L) {
			if(is.name(substitute(global.model))) {
				cry(, "call stored in 'global.model' contains dotted names and cannot be updated. \n    Consider using 'updateable' on the modelling function")
			} else gmCall[isDotted] <-
				substitute(global.model)[names(gmCall[isDotted])]
		}
		# object from 'run.mark.model' has $call of 'make.mark.model' - fixing
		# it here:
		if(inherits(global.model, "mark") && gmCall[[1L]] == "make.mark.model") {
			gmCall <- call("run.mark.model", model = gmCall, invisible = TRUE)
		}
	}

	lik <- .getLik(global.model)
	logLik <- lik$logLik

	# *** Rank ***
	rank.custom <- !missing(rank)

	if(!rank.custom && lik$name == "qLik") {
		rank <- "QIC"
		cry(, "using 'QIC' instead of 'AICc'", warn = TRUE)
	}

	rankArgs <- list(...)

	if(any(badargs <- names(rankArgs) == "marg.ex")) {
		cry(, "argument \"marg.ex\" is defunct and has been ignored",
			 warn = TRUE)
		rankArgs <- rankArgs[!badargs]
	}
	if(any(names(rankArgs) == "na.action"))
		cry("RTFM", "argument \"na.action\" is inappropriate here",
			 warn = FALSE)

	IC <- .getRank(rank, rankArgs)

	if(any(badargs <- is.na(match(names(rankArgs),
		c(names(formals(get("rank", environment(IC))))[-1L], names(formals()))))))
		cry("RTFM", ngettext(sum(badargs),
			"argument %s is not a name of formal argument of %s",
			"arguments %s are not names of formal arguments of %s"),
			prettyEnumStr(names(rankArgs[badargs])), "'pdredge' or 'rank'",
			warn = TRUE)

	ICName <- as.character(attr(IC, "call")[[1L]])

	if(length(tryCatch(IC(global.model), error = function(e) {
		stop(simpleError(conditionMessage(e), subst(attr(IC, "call"),
			x = as.name("global.model"))))
	})) != 1L) {
		cry(, "result of '%s' is not of length 1", asChar(attr(IC, "call")))
	}

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nIntercepts <- sum(attr(allTerms, "intercept"))


###PAR
	# parallel: check whether the models would be identical:
	if(doParallel && check) testUpdatedObj(cluster, global.model, gmCall, level = check)
###PAR

	# Check for na.omit
	if(!(gmNaAction <- .checkNaAction(cl = gmCall, what = "'global.model'")))
		cry(, attr(gmNaAction, "message"))


	if(names(gmCall)[2L] == "") gmCall <-
		match.call(gmCall, definition = eval.parent(gmCall[[1L]]),
				   expand.dots = TRUE)

    gmCoefNames <- names(coeffs(global.model))
    if(any(dup <- duplicated(gmCoefNames)))
        cry(, "model cannot have duplicated coefficient names: ",
             prettyEnumStr(gmCoefNames[dup]))

	gmCoefNames <- fixCoefNames(gmCoefNames)

	nVars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		cry(, "comparing models fitted by REML", warn = TRUE)

	if ((betaMode != 0L) && is.null(tryCatch(std.coef(global.model, betaMode == 2L),
		error = return_null, warning = return_null))) {
		cry(, "do not know how to standardize coefficients of '%s', argument 'beta' ignored",
			 class(global.model)[1L], warn = TRUE)
		betaMode <- 0L
		strbeta <- "none"
	}

	if(nomlim <- is.null(m.lim)) m.lim <- c(0, NA)
	## XXX: backward compatibility:
	if(!missing(m.max) || !missing(m.min)) {
		warning("arguments 'm.min' and 'm.max' are deprecated, use 'm.lim' instead")
		if(!nomlim) stop("cannot use both 'm.lim' and 'm.min' or 'm.max'")
		if(!missing(m.min)) m.lim[1L] <- m.min[1L]
		if(!missing(m.max)) m.lim[2L] <- m.max[1L]
	}
	if(!is.numeric(m.lim) || length(m.lim) != 2L || any(m.lim < 0, na.rm = TRUE))
		stop("invalid 'm.lim' value")
	m.lim[2L] <- if (!is.finite(m.lim[2L])) (nVars - nIntercepts) else
		min(nVars - nIntercepts, m.lim[2L])
	if (!is.finite(m.lim[1L])) m.lim[1L] <- 0
	m.min <- m.lim[1L]
    m.max <- m.lim[2L]

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				cry(, "'fixed' should be a one-sided formula", warn = TRUE)
			fixed <- as.vector(getAllTerms(fixed))
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
		} else if (!is.character(fixed)) {
			cry(, paste("'fixed' should be either a character vector with",
						   " names of variables or a one-sided formula"))
		}
		if (!all(i <- (fixed %in% allTerms))) {
			cry(, "some terms in 'fixed' do not exist in 'global.model': %s",
				 prettyEnumStr(fixed[!i]), warn = TRUE)
			fixed <- fixed[i]
		}
	}

	deps <- attr(allTerms0, "deps")
	fixed <- union(fixed, rownames(deps)[rowSums(deps, na.rm = TRUE) == ncol(deps)])
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])

	nFixed <- length(fixed)
	if(nFixed != 0L) message(sprintf(ngettext(nFixed, "Fixed term is %s", "Fixed terms are %s"),
		prettyEnumStr(fixed)))

	termsOrder <- order(allTerms %in% fixed)
	allTerms <- allTerms[termsOrder]

	di <- match(allTerms, rownames(deps))
	deps <- deps[di, di]

	gmFormulaEnv <- environment(as.formula(formula(global.model), env = gmEnv))
	# TODO: gmEnv <- gmFormulaEnv ???

	### BEGIN Manage 'varying'
	## @param:	varying
	## @value:	varying, varyingNames, variants, nVariants, nVarying
	if(!missing(varying) && !is.null(varying)) {
		nVarying <- length(varying)
		varyingNames <- names(varying)
		fvarying <- unlist(varying, recursive = FALSE, use.names = FALSE)
		vlen <- vapply(varying, length, 1L)
		nVariants <- prod(vlen)
		variants <- as.matrix(expand.grid(split(seq_len(sum(vlen)),
			rep(seq_len(nVarying), vlen))))

		variantsFlat <- unlist(lapply(varying, .makeListNames),
			recursive = FALSE, use.names = FALSE)

	} else {
		variants <- varyingNames <- NULL
		nVariants <- 1L
		nVarying <- 0L
	}
	## END: varying

	## BEGIN Manage 'extra'
	## @param:	extra, global.model, gmFormulaEnv,
	## @value:	extra, nextra, extraNames, nullfit_
	if(!missing(extra) && length(extra) != 0L) {
		# a cumbersome way of evaluating a non-exported function in a parent frame:
		extra <- eval(as.call(list(call("get", ".get.extras",
			envir = call("asNamespace", .packageName), inherits = FALSE),
					 substitute(extra), r2nullfit = TRUE)), parent.frame())

		#extra <- eval(call(".get.extras", substitute(extra), r2nullfit = TRUE), parent.frame())
		if(any(c("adjR^2", "R^2") %in% names(extra))) {
			nullfit_ <- null.fit(global.model, evaluate = TRUE, envir = gmFormulaEnv)
		}
		applyExtras <- function(x) unlist(lapply(extra, function(f) f(x)))
		extraResult <- applyExtras(global.model)
		if(!is.numeric(extraResult))
			cry(, "function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## END: manage 'extra'

	nov <- as.integer(nVars - nFixed)
	ncomb <- (2L ^ nov) * nVariants

	if(nov > 31L) cry(, "number of predictors [%d] exceeds allowed maximum of 31", nov)
	#if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f combinations", nov, ncomb))
	nmax <- ncomb * nVariants
	rvChunk <- 25L
	if(evaluate) {
		rvNcol <- nVars + nVarying + 3L + nextra
		rval <- matrix(NA_real_, ncol = rvNcol, nrow = rvChunk)
		coefTables <- vector(rvChunk, mode = "list")
	}


	## BEGIN: Manage 'subset'
	## @param:	hasSubset, subset, allTerms, [interceptLabel],
	## @value:	hasSubset, subset
	if(missing(subset))  {
		hasSubset <- 1L
	} else {
		if(!tryCatch(is.language(subset) || is.matrix(subset), error = function(e) FALSE))
			subset <- substitute(subset)

		if(is.matrix(subset)) {
			dn <- dimnames(subset)
			#at <- allTerms[!(allTerms %in% interceptLabel)]
			n <- length(allTerms)
			if(is.null(dn) || any(sapply(dn, is.null))) {
				di <- dim(subset)
				if(any(di != n)) stop("unnamed 'subset' matrix does not have both dimensions",
					" equal to number of terms in 'global.model': %d", n)

				dimnames(subset) <- list(allTerms, allTerms)
			} else {
				if(!all(unique(unlist(dn)) %in% allTerms))
					warning("at least some dimnames of 'subset' matrix do not ",
					"match term names in 'global.model'")

				subset0 <- subset
				subset <- matrix(subset[
					match(allTerms, rownames(subset)),
					match(allTerms, colnames(subset))],
					dimnames = list(allTerms, allTerms),
					nrow = n, ncol = n)
				nas <- is.na(subset)
				lotri <- lower.tri(subset)
				i <- lotri & nas & !t(nas)
				subset[i] <- t(subset)[i]
				subset[!lotri] <- NA

			}
			if(any(!is.na(subset[!lower.tri(subset)]))) {
				warning("non-missing values exist outside the lower triangle of 'subset'")
				subset[!lower.tri(subset)] <- NA
			}
			mode(subset) <- "logical"
			hasSubset <- 2L # subset as matrix
		} else {
			if(inherits(subset, "formula")) {
				if (subset[[1L]] != "~" || length(subset) != 2L)
					stop("'subset' formula should be one-sided")
				subset <- subset[[2L]]
			}
			subset <- as.expression(subset)
			ssValidNames <- c("comb", "*nvar*")

			tmpTerms <- terms(reformulate(allTerms0[!(allTerms0 %in% interceptLabel)]))
			gloFactorTable <- t(attr(tmpTerms, "factors") != 0)

			offsetNames <- sapply(attr(tmpTerms, "variables")[attr(tmpTerms, "offset") + 1L], asChar)
			if(length(offsetNames) != 0L) {
				gloFactorTable <- rbind(gloFactorTable,
					matrix(FALSE, ncol = ncol(gloFactorTable), nrow = length(offsetNames),
						dimnames = list(offsetNames, NULL)))
				for(i in offsetNames) gloFactorTable[offsetNames, offsetNames] <- TRUE
				#Note `diag<-` does not work for x[1x1] matrix:
				# diag(gloFactorTable[offsetNames, offsetNames, drop = FALSE]) <- TRUE
			}

			DebugPrint(gloFactorTable)

			# fix interaction names in rownames:
			rownames(gloFactorTable) <- allTerms0[!(allTerms0 %in% interceptLabel)]

			subsetExpr <- subset[[1L]]
			subsetExpr <- exprapply0(subsetExpr, ".", .sub_dot, gloFactorTable,
				allTerms, as.name("comb"))
			subsetExpr <- exprapply0(subsetExpr, c("{", "Term"), .sub_Term)

			tmp <- updateDeps(subsetExpr, deps)
			subsetExpr <- tmp$expr
			deps <- tmp$deps

			subsetExpr <- exprapply0(subsetExpr, "dc", .sub_args_as_vars)
			subsetExpr <- .subst4Vec(subsetExpr, allTerms, "comb")

			if(nVarying) {
				ssValidNames <- c("cVar", "comb", "*nvar*")
					subsetExpr <- exprapply0(subsetExpr, "V", .sub_V,
						as.name("cVar"), varyingNames)
				if(!all(all.vars(subsetExpr) %in% ssValidNames))
					subsetExpr <- .subst4Vec(subsetExpr, varyingNames,
						"cVar", fun = "[[")
			}
			ssVars <- all.vars(subsetExpr)
			okVars <- ssVars %in% ssValidNames
			if(!all(okVars)) stop("unrecognized names in 'subset' expression: ",
				prettyEnumStr(ssVars[!okVars]))

			ssEnv <- new.env(parent = parent.frame())
			ssFunc <- setdiff(all.vars(subsetExpr, functions = TRUE), ssVars)
			if("dc" %in% ssFunc) assign("dc", .subset_dc, ssEnv)

			hasSubset <- if(any(ssVars == "cVar")) 4L else # subset as expression
				3L # subset as expression using 'varying' variables

		}
	} # END: manage 'subset'

	comb.sfx <- rep(TRUE, nFixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	extraResult1 <- integer(0L)
	calls <- vector(mode = "list", length = rvChunk)
	ord <- integer(rvChunk)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nIntercepts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		gmDataHead = if(!is.null(gmCall$data)) {
			if(eval(call("is.data.frame", gmCall$data), gmEnv))
				eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
			} else NULL,
		gmFormulaEnv = gmFormulaEnv
		)


	# BEGIN parallel
	qi <- 0L
	queued <- vector(qlen, mode = "list")
	props <- list(
				gmEnv = gmEnv,
				IC = IC,
				# beta = beta,
				# allTerms = allTerms,
				nextra = nextra,
				matchCoefCall = as.call(c(list(
					as.name("matchCoef"), as.name("fit1"),
					all.terms = allTerms, beta = betaMode,
					allCoef = TRUE), ct.args))
				# matchCoefCall = as.call(c(alist(matchCoef, fit1, all.terms = Z$allTerms,
				#   beta = Z$beta, allCoef = TRUE), ct.args))
		)
	if(nextra) {
		props$applyExtras <- applyExtras
		props$extraResultNames <- names(extraResult)
	}
	props <- as.environment(props)

	if(doParallel) {
		clusterVExport(cluster,   pdredge_props = props,
								  .pdredge_process_model = pdredge_process_model
								  )
		clusterCall(cluster, eval, call("options", options("na.action")), env = 0L)
	}
	# END parallel

	retColIdx <- if(nVarying) -nVars - seq_len(nVarying) else TRUE

	if(trace > 1L) {
		progressBar <- if(.Platform$GUI == "Rgui") {
			 utils::winProgressBar(max = ncomb, title = "'dredge' in progress")
		} else utils::txtProgressBar(max = ncomb, style = 3)
		setProgressBar <- switch(class(progressBar),
			   txtProgressBar = utils::setTxtProgressBar,
			   winProgressBar = utils::setWinProgressBar,
			function(...) {})
		on.exit(close(progressBar))
	}


	warningList <- list()

	iComb <- -1L
	while((iComb <- iComb + 1L) < ncomb) {
		varComb <- iComb %% nVariants
		jComb <- (iComb - varComb) / nVariants
		#print(c(iComb, jComb, ncomb, varComb + 1L))
		if(varComb == 0L) {
			isok <- TRUE

			comb <- c(as.logical(intToBits(jComb)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nIntercepts

			# !!! POSITIVE condition for 'pdredge', NEGATIVE for 'dredge':
			if((nvar >= m.min && nvar <= m.max) &&
				formula_margin_check(comb, deps) &&
				switch(hasSubset,
					# 1 - no subset, 2 - matrix, 3 - expression
					TRUE,                                    # 1
					all(subset[comb, comb], na.rm = TRUE),   # 2
					evalExprInEnv(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),		 # 3
					TRUE
					)
				) {

				newArgs <- makeArgs(global.model, allTerms[comb], argsOptions) #comb
				formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
					attr(newArgs, "formulaList")

				if(!is.null(attr(newArgs, "problems"))) {
					print.warnings(structure(vector(mode = "list",
						length = length(attr(newArgs, "problems"))),
							names = attr(newArgs, "problems")))
				} # end if <problems>
				cl <- gmCall
				cl[names(newArgs)] <- newArgs
			} else isok <- FALSE # end if <subset, m.max >= nvar >= m.min>
		} #  end if(jComb != prevJComb)

		if(isok) {
			## --- Variants ---------------------------
			clVariant <- cl
			isok2 <- TRUE
			if(nVarying) {
				cvi <- variants[varComb + 1L, ]
				isok2 <- (hasSubset != 4L) || evalExprInEnv(subsetExpr, env = ssEnv,
					enclos = parent.frame(), comb = comb, `*nvar*` = nvar,
					cVar = variantsFlat[cvi])
				clVariant[varyingNames] <- fvarying[cvi]
			}

			if(isok2) {
				if(evaluate) {
					if(trace == 1L) {
						cat(iComb, ": "); print(clVariant)
						utils::flush.console()
					} else if(trace == 2L) {
						setProgressBar(progressBar, value = iComb,
							title = sprintf("pdredge: %d of %.0f subsets", k, (k / iComb) * ncomb))
					}
					qi <- qi + 1L
					queued[[(qi)]] <- list(call = clVariant, id = iComb)
				} else { # if !evaluate
					k <- k + 1L # all OK, add model to table
					rvlen <- length(ord)
					if(k > rvlen) {
						nadd <- min(rvChunk, nmax - rvlen)
						#message(sprintf("extending result from %d to %d", rvlen, rvlen + nadd))
						addi <- seq.int(rvlen + 1L, length.out = nadd)
						calls[addi] <- vector("list", nadd)
						ord[addi] <- integer(nadd)
					}
					calls[[k]] <- clVariant
					ord[k] <- iComb
				}
			}
		} # if isok

		#if(evaluate && qi && (qi + nvariants > qlen || iComb == ncomb)) {
		if(evaluate && qi && (qi > qlen || (iComb + 1L) == ncomb)) {
			qseq <- seq_len(qi)
			qresult <- .getRow(queued[qseq])
			utils::flush.console()

			if(any(vapply(qresult, is.null, TRUE)))
				stop("some results returned from cluster node(s) are NULL. \n",
					"This should not happen and indicates problems with ",
					"the cluster node", domain = "R-MuMIn")
			haveProblems <- logical(qi)

			nadd <- sum(sapply(qresult, function(x) inherits(x$value, "condition")
				+ length(x$warnings)))
			wi <- length(warningList)
			if(nadd) warningList <- c(warningList, vector(nadd, mode = "list"))

			# DEBUG: print(sprintf("Added %d warnings, now is %d", nadd, length(warningList)))

			for (i in qseq)
				for(cond in c(qresult[[i]]$warnings,
					if(inherits(qresult[[i]]$value, "condition"))
						list(qresult[[i]]$value))) {
						wi <- wi + 1L
						warningList[[wi]] <- if(is.null(conditionCall(cond)))
							queued[[i]]$call else conditionCall(cond)
						if(inherits(cond, "error")) {
							haveProblems[i] <- TRUE
							msgsfx <- "(model %d skipped)"
						} else
							msgsfx <- "(in model %d)"
						names(warningList)[wi] <- paste(conditionMessage(cond),
							 gettextf(msgsfx, queued[[i]]$id))
						attr(warningList[[wi]], "id") <- queued[[i]]$id
				}

			withoutProblems <- which(!haveProblems)
			qrows <- lapply(qresult[withoutProblems], "[[", "value")
			qresultLen <- length(qrows)
			rvlen <- nrow(rval)

			if(retNeedsExtending <- k + qresultLen > rvlen) {
				nadd <- min(max(rvChunk, qresultLen), nmax - rvlen)
				rval <- rbind(rval, matrix(NA_real_, ncol = rvNcol, nrow = nadd),
					deparse.level = 0L)
				addi <- seq.int(rvlen + 1L, length.out = nadd)
				coefTables[addi] <- vector("list", nadd)
				calls[addi] <- vector("list", nadd)
				ord[addi] <- integer(nadd)
			}
			qseqOK <- seq_len(qresultLen)
			for(m in qseqOK) rval[k + m, retColIdx] <- qrows[[m]]
			ord[k + qseqOK] <- vapply(queued[withoutProblems], "[[", 1L, "id")
			calls[k + qseqOK] <- lapply(queued[withoutProblems], "[[", "call")
			coefTables[k + qseqOK] <- lapply(qresult[withoutProblems], "[[", "coefTable")
			k <- k + qresultLen
			qi <- 0L
		}
	} ### for (iComb ...)

	if(k == 0L) {
		if(length(warningList)) print.warnings(warningList)
		stop("the result is empty")
	}

	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(rval)) {
		i <- seq_len(k)
		rval <- rval[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nVarying) {
		varlev <- ord %% nVariants
		varlev[varlev == 0L] <- nVariants
		rval[, nVars + seq_len(nVarying)] <- variants[varlev, ]
	}

	rval <- as.data.frame(rval)
	row.names(rval) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))
	rval[tfac] <- lapply(rval[tfac], factor, levels = NaN, labels = "+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	rval[, i] <- rval[, v]
	allTerms <- allTerms[v]
	colnames(rval) <- c(allTerms, varyingNames, extraNames, "df", lik$name, ICName)

	if(nVarying) {
		variant.names <- vapply(variantsFlat, asChar, "", width.cutoff = 20L)
		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nVarying), vlen))
		names(vnum) <- varyingNames
		for (i in varyingNames) rval[, i] <-
			factor(rval[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])

	}

	rval <- rval[o <- order(rval[, ICName], decreasing = FALSE), ]
	coefTables <- coefTables[o]

	rval$delta <- rval[, ICName] - min(rval[, ICName])
	rval$weight <- exp(-rval$delta / 2) / sum(exp(-rval$delta / 2))
    mode(rval$df) <- "integer"

	rval <- 
    structure(rval,
		model.calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		beta = strbeta,
		call = match.call(expand.dots = TRUE),
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varyingNames, ## XXX: remove
        column.types = {
			colTypes <- c(terms = length(allTerms), varying = length(varyingNames),
				extra = length(extraNames), df = 1L, loglik = 1L, ic = 1L, delta = 1L,
				weight = 1L)
			column.types <- rep(1L:length(colTypes), colTypes)
			names(column.types) <- colnames(rval)
			lv <- 1L:length(colTypes)
			factor(column.types, levels = lv, labels = names(colTypes)[lv])
		},
        class = c("model.selection", "data.frame")
	)

	if(length(warningList)) {
		class(warningList) <- c("warnings", "list")
		attr(rval, "warnings") <- warningList
	}

	if (!is.null(attr(allTerms0, "random.terms")))
		attr(rval, "random.terms") <- attr(allTerms0, "random.terms")

	if(doParallel) clusterCall(cluster, "rm",
		list = c(".pdredge_process_model", "pdredge_props"), envir = .GlobalEnv)


	return(rval)
} ######



`pdredge_process_model` <- function(modv, envir = get("pdredge_props", .GlobalEnv)) {
	### modv == list(call = clVariant, id = modelId)
	result <- tryCatchWE(eval(modv$call, get("gmEnv", envir)))
	if (inherits(result$value, "condition")) return(result)

	fit1 <- result$value
	if(get("nextra", envir) != 0L) {
		extraResult1 <- get("applyExtras", envir)(fit1)
		nextra <- get("nextra", envir)
		if(length(extraResult1) < nextra) {
			tmp <- rep(NA_real_, nextra)
			tmp[match(names(extraResult1), get("extraResultNames", envir))] <-
				extraResult1
			extraResult1 <- tmp
		}
	} else extraResult1 <- NULL
	ll <- .getLik(fit1)$logLik(fit1)

	#mcoef <- matchCoef(fit1, all.terms = get("allTerms", envir),
	# beta = get("beta", envir), allCoef = TRUE)
	mcoef <- eval(get("matchCoefCall", envir))

	list(value = c(mcoef, extraResult1, df = attr(ll, "df"), ll = ll,
		ic = get("IC", envir)(fit1)),
		nobs = nobs(fit1),
		coefTable = attr(mcoef, "coefTable"),
		warnings = result$warnings)

}

.test_pdredge <- function(dd) {
	cl <- attr(dd, "call")
	cl$cluster <- cl$check <- NULL
	cl[[1L]] <- as.name("dredge")
	if(!identical(c(dd), c(eval(cl)))) stop("Whoops...")
	dd
}
