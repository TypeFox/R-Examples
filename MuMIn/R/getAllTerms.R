`getAllTerms.default` <-
#function(x, ...) getAllTerms.formula(as.formula(formula(x)), ...)
function(x, ...) getAllTerms.terms(terms(as.formula(formula(x))), ...)

`getAllTerms.gam` <-
function(x, intercept = FALSE, ...)
	getAllTerms.terms(terms(formula(x), ...), intercept = intercept)

`getAllTerms.lm` <-
function(x, intercept = FALSE, ...)
	getAllTerms.terms(terms(x, ...), intercept = intercept)

`getAllTerms.terms` <-
function(x, offset = TRUE, intercept = FALSE, ...) {

	interceptLabel <- "(Intercept)"
	variables <- attr(x, "variables")[-1L]

	if (!is.null(attr(x, "offset"))){
		offs <- sapply(variables[attr(x, "offset")], deparse)
	} else offs <- NULL

	ans <- attr(x, "term.labels")

	# Get term names, with higher order term components arranged alphabetically
	if (length(ans) > 0L) {
		factors <- attr(x, "factors")
		factors <- factors[order(rownames(factors)), , drop = FALSE]
		v <- rownames(factors)
		ans <- apply(factors != 0L, 2L, function(x) paste0(v[x], collapse = ":"))
	}

	# Leave out random terms (lmer type)
	#ran <- attr(x, "variables")[-1][-c(attr(x, "offset"), attr(x, "response"))]
	ran <- as.character(variables[vapply(variables, function(x) length(x) == 3L && x[[1L]] == "|", TRUE)])
	ifx <- !(ans %in% ran)

	ans <- ans[ifx] # ifx - indexes of fixed terms
	#retUnsorted <- ans

	# finally, sort by term order and then alphabetically
	#ans <- unname(ans[order(attr(x, "order")[ifx], ans)])
	ord <- order(attr(x, "order")[ifx], gsub("I\\((.*)\\)", "\\1", ans))
	ans <- unname(ans[ord])

	deps <- if (length(ans) > 0L) termdepmat(reformulate(ans)) else
		matrix(FALSE, 0L, 0L)
		
	dimnames(deps) <- list(ans, ans)
	diag(deps) <- NA
	
	if(intercept && attr(x, "intercept")) {
		ans <- c(interceptLabel, ans)
		ord <- c(1L, ord + 1L)
	}

	if (!is.null(offs[1L])) {
		if (offset) {
			ans <- c(ans, offs)
			ord <- c(ord, length(ord) + 1L)
		}
		attr(ans, "offset") <- offs
	}
	attr(ans, "intercept") <- attr(x, "intercept")
	attr(ans, "interceptLabel") <- interceptLabel

	if (length(ran) > 0L) {
		attr(ans, "random.terms") <- ran
		f.random <- reformulate(c(".", paste0("(", ran, ")")), response = ".")
		environment(f.random) <- environment(x)
		attr(ans, "random") <- f.random
	}

	response <- attr(x, "response")
	response <- if(response == 0L) NULL else variables[[response]]
	attr(ans, "response") <- response
	attr(ans, "order") <- order(ord)
	attr(ans, "deps") <- deps
	ans
}

`getAllTerms.formula` <-
function(x, ...) getAllTerms.terms(terms.formula(x), ...)

`getAllTerms.lme` <-
function(x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	attr(ret, "random") <- . ~ .

	# Code from nlme:::print.reStruct, modified slightly
	reStruct <- x$modelStruct$reStruct
	nobj <- length(reStruct)
	if (is.null(namx <- names(reStruct)))
		names(reStruct) <- nobj:1L
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	reStruct[] <- rev(reStruct)
	aux <- t(array(rep(names(reStruct), nobj), c(nobj, nobj)))
	aux[lower.tri(aux)] <- ""
	attr(ret, "random.terms") <- paste(lapply(lapply(reStruct, attr, "formula"),
		"[[", 2L), "|",
		rev(apply(aux, 1L, function(z) paste(z[z != ""], collapse = " %in% "))))

	return(ret)
}

# Apparently there is no (explicit) intercept in coxph, but 'terms' gives
# attr(,"intercept") == 1.
`getAllTerms.coxph` <- function (x, ...) {
	ret <- getAllTerms.default(x, ...)
	attr(ret, "intercept") <- 0L
	attr(ret, "interceptLabel") <- NULL
	return(ret)
}

`getAllTerms.glmmML` <- function (x, ...) {
	ret <- getAllTerms.terms(terms(x), ...)
	attr(ret, "random.terms") <-  paste("1 |",  x$call$cluster)
	return(ret)
}

#`getAllTerms.hurdle` <- function(x, intercept = FALSE, ...) {
#	f <- as.formula(formula(x))
#	# to deal with a dot in formula (other classes seem to expand it)
#	if("." %in% all.vars(f))
#		getAllTerms.terms(terms.formula(f, data = eval(x$call$data, envir = environment(f)))
#			
#			, intercept = intercept)
#	else getAllTerms.formula(f, intercept = intercept)
#}

split_formula_by_bar <- function(f) {
	n <- length(f)
	ans <- if(length(f[[n]]) != 1L && f[[n]][[1L]] == "|") {
		f1 <- vector("list", 2L)
		for(i in 1L:2L) {
			f1[[i]] <- f
			f1[[i]][[n]] <- f[[n]][[i+ 1]]
		}
		f1
	} else list(f)
	ans
}


`getAllTerms.hurdle` <- 
`getAllTerms.zeroinfl` <-
function(x, intercept = FALSE, ...) {

	formList <- split_formula_by_bar(formula(x))
	formList <- lapply(lapply(formList, terms.formula, data = eval(x$call$data)),
		formula)
	z <- lapply(formList, getAllTerms, intercept = TRUE)
	
	if(oneform <- length(formList) == 1L) z <- c(z, z)

	deps <- termdepmat_combine(lapply(z, attr, "deps"))

	ord <- unlist(lapply(z, attr, "order"))
	n <- sapply(z, length)
	if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
		
	zz <- unlist(z)
	interceptIdx <- zz == "(Intercept)"
	offsetIdx <- match(zz, unique(unlist(lapply(z, attr,"offset"))), nomatch = 0) != 0
	termIdx <- !(offsetIdx | interceptIdx)

	zz <- paste0(rep(c("count", "zero")[seq_along(z)], sapply(z, length)),
				 "_", zz)
	
	dimnames(deps) <- list(zz[termIdx], zz[termIdx])
	
	if(oneform) { # dependency of count_X and zero_X
		k <- length(zz[termIdx]) / 2
		deps[c(seq(k + 1L, by = 2L * k + 1L, length.out = k),
			   seq((2L * k * k) + 1L, by = 2L * k + 1L, length.out = k))] <- TRUE
	}
	
	ret <- if(!intercept) zz[!interceptIdx] else zz
	if(any(offsetIdx)) attr(ret, "offset") <- zz[offsetIdx]
	attr(ret, "intercept") <- pmin(which(interceptIdx), 1)
	attr(ret, "interceptLabel") <- zz[interceptIdx]
	attr(ret, "response") <- attr(z[[1L]], "response")
	attr(ret, "order") <- if(!intercept) order(ord[!interceptIdx]) else ord
	attr(ret, "deps") <- deps
	ret
}

## TODO: test with offsets
`getAllTerms.betareg` <-
function(x, intercept = FALSE, ...) {
	formList <- split_formula_by_bar(formula(x))
	formList <- lapply(lapply(formList, terms.formula, data = model.frame(x)),
		formula)
	oneform <- length(formList) == 1L
	z <- lapply(formList, getAllTerms, intercept = TRUE)
	
	deps <- termdepmat_combine(lapply(z, attr, "deps"))
	
	ord <- unlist(lapply(z, attr, "order"))
	n <- sapply(z, length)
	if(length(z) > 1L) ord[-j] <- ord[-(j <- seq_len(n[1L]))] + n[1L]
	zz <- unlist(z)
	interceptIdx <- zz == "(Intercept)"
	offsetIdx <- match(zz, unique(unlist(lapply(z, attr,"offset"))), nomatch = 0) != 0
	termIdx <- !(offsetIdx | interceptIdx)
	
	if(!oneform && n[2L] != 0L) {
		i.phi <- -seq.int(n[1L])
		zz[i.phi] <- paste("(phi)", zz[i.phi], sep = "_")
	}
	dimnames(deps) <- list(zz[termIdx], zz[termIdx])
	
	ret <- if(!intercept) zz[!interceptIdx] else zz
	if(any(offsetIdx)) attr(ret, "offset") <- zz[offsetIdx]
	attr(ret, "intercept") <- pmin(which(interceptIdx), 1)
	attr(ret, "interceptLabel") <- zz[interceptIdx]
	attr(ret, "response") <- attr(z[[1L]], "response")
	attr(ret, "order") <- if(!intercept) order(ord[!interceptIdx]) else ord
	attr(ret, "deps") <- deps
	ret
}




`getAllTerms.glimML` <- function(x, intercept = FALSE, ...) {
	ret <- getAllTerms.default(x, intercept = intercept, ...)
	ttran <- terms.formula(x@random)
	ran <- attr(ttran, "term.labels")
	if(length(ran)) attr(ret, "random.terms") <- paste("1 |", ran)
	ret
}

`getAllTerms.coxme` <-
function(x, ...)  {
	ret <- getAllTerms.terms(terms(x))
	random <- x$formulaList$random
	attr(ret, "random.terms") <- as.character(random)
	f <- as.name(".")
	for(f1 in random) f <- call("+", f, f1)
	attr(ret, "random") <- call("~", as.name("."), f)
	attr(ret, "intercept") <- 0L
	attr(ret, "interceptLabel") <- NULL
	ret
}


#gdistsamp -> unmarkedFitGDS
#lambda = abundance / lambda
#phi = availability / alpha
#p = detection / det

get_all_terms_multiple_form <-
function(f, fnames, intercept, ...) {
	ret <- vector("list", length(fnames))
	i <- 0L
	while(is.call(f) && f[[1L]] == "~") {
		ret[[i <- i + 1L]] <- as.formula(f[c(1L, length(f))])
		f <- f[[2L]]
	}
	ret <- lapply(rev(ret), `environment<-`, NULL)
	names(ret) <- fnames
	ret <- lapply(ret, getAllTerms.formula, intercept = FALSE)
	
	deps <- termdepmat_combine(lapply(ret, attr, "deps"))

	attrInt <- sapply(ret, attr, "intercept")
	ret <- unlist(lapply(names(ret), function(i) if(length(ret[[i]]))
						 paste0(i, "(", ret[[i]], ")") else character(0L)))

	dimnames(deps) <- list(ret, ret)

	Ints <- paste0(names(attrInt[attrInt != 0L]), "(Int)")
	if(intercept) ret <- c(Ints, ret)
	attr(ret, "intercept") <- attrInt
	attr(ret, "interceptLabel") <- Ints
	attr(ret, "deps") <- deps
	return(ret)
}

`getAllTerms.unmarkedFitGDS` <- function (x, intercept = FALSE, ...)  {
	get_all_terms_multiple_form(formula(x), c("lambda", "alpha", "det"), intercept = intercept, ...)
}

## TODO:
#
#`getAllTerms.unmarkedFit` <- function (x, intercept = FALSE, ...)  {
#	get_all_terms_multiple_form(formula(x), ......)
#}


`getAllTerms.unmarkedFit` <- function (x, intercept = FALSE, ...)  {
	f <- formula(x)
	ret <- list()
	while(is.call(f) && f[[1L]] == "~") {
		ret <- c(ret, as.formula(f[c(1L, length(f))]))
		f <- f[[2L]]
	}
	ret <- lapply(ret, `environment<-`, NULL)
	names(ret) <- sapply(x@estimates@estimates, slot, "short.name")[seq_along(ret)]
	
	ret <- lapply(ret, getAllTerms.formula, intercept = FALSE)
	
	deps <- termdepmat_combine(lapply(ret, attr, "deps"))

	attrInt <- sapply(ret, attr, "intercept")
	#ret <- unlist(lapply(names(ret), function(i) sprintf("%s(%s)", i, ret[[i]])))
	ret <- unlist(lapply(names(ret), function(i) if(length(ret[[i]]))
						 paste0(i, "(", ret[[i]], ")") else character(0L)))

	dimnames(deps) <- list(ret, ret)

	Ints <- paste0(names(attrInt[attrInt != 0L]), "(Int)")
	if(intercept) ret <- c(Ints, ret)
	attr(ret, "intercept") <- attrInt
	attr(ret, "interceptLabel") <- Ints
	attr(ret, "deps") <- deps
	return(ret)
}

## tweak for 'distsamp' models: prefix the detection "p(...)" terms with 'sigma'
`getAllTerms.unmarkedFitDS` <- function (x, intercept = FALSE, ...)  {
	tt <- getAllTerms.unmarkedFit(x, intercept = FALSE)
	ret <- gsub("^p\\(", "p(sigma", c(tt))
	intLab <- attr(tt, "interceptLabel")
	intLab[intLab == "p(Int)"] <- "p(sigma(Intercept))"
	if(intercept) ret <- c(intLab, ret)
	mostattributes(ret) <- attributes(tt)
	attr(ret, "interceptLabel") <- intLab
	deps <- attr(ret, "deps")
	dn <- gsub("^p\\(", "p(sigma", rownames(deps))
	dimnames(deps) <- list(dn, dn)
	attr(ret, "deps") <- deps
	ret
}

`getAllTerms.MCMCglmm` <- 
function (x, ...) {
	res <- getAllTerms.default(x, ...) 
	attr(res, "random") <- .formulaEnv(.~., environment(formula(x)))
	attr(res, "random.terms") <- asChar(x$Random$formula)[1L]
	res
}

`getAllTerms.gamm` <-
function (x, ...) getAllTerms(x$gam, ...)


`getAllTerms.mark` <- 
function (x, intercept = FALSE, ...) {
	
	f <- formula(x, expand = FALSE)[[2L]]
	formlist <- list()
	while(length(f) == 3L && f[[1L]] == "+") {
		formlist <- c(f[[3L]], formlist)
		f <- f[[2L]]
	}
	formlist <- append(f, formlist)
	
	wrapfunc <- function(x, func) if(length(x) == 0L) x else paste0(func, "(", x, ")")

	alltermlist <- lapply(formlist, function(x, intercept) {
		func <- asChar(x[[1L]])
		at <- getAllTerms(terms(eval(call("~", x[[2L]]))), intercept = intercept)
		at[] <- wrapfunc(at, func)
		dn <- wrapfunc(rownames(attr(at, "deps")), func)
		attr(at, "interceptLabel") <- wrapfunc(attr(at, "interceptLabel"), func)
		dimnames(attr(at, "deps")) <- list(dn, dn)
		at
	}, intercept)
	
	retval <- unlist(alltermlist, recursive = TRUE)
	for(a in c("intercept", "interceptLabel")) {
		attr(retval, a) <-	unlist(sapply(alltermlist, attr, a))
	}
	attr(retval, "order") <- order(rep(seq_along(alltermlist), vapply(alltermlist, length, 1L)),
		unlist(lapply(alltermlist, attr, "order")))
	attr(retval, "deps") <- termdepmat_combine(lapply(alltermlist, attr, "deps"))
	retval
}

`getAllTerms.asreml`  <-
function(x, intercept = FALSE, ...)
getAllTerms.terms(terms(formula(x), ...), intercept = intercept)

`getAllTerms.cpglmm` <-
function (x, intercept = FALSE, ...) 
getAllTerms(x@formula, intercept = intercept)


`getAllTerms` <-
function(x, ...)
UseMethod("getAllTerms")

# TODO: return object of class 'allTerms'
print.allTerms <-
function(x, ...) {
	cat("Model terms: \n")
	if(!length(x)) {
		cat("<None> \n")
	} else {
		print.default(as.vector(x), quote = TRUE)
	}
	ints <- attr(x, "interceptLabel")
	if(!is.null(ints)) {
		cat(ngettext(n = length(ints), "Intercept:", "Intercepts:"), "\n")
		print.default(ints,quote = TRUE)
	}
}
