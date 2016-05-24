`null.fit` <-
function(x, evaluate = FALSE, RE.keep = FALSE, envir = NULL) {
	cl <- get_call(x)
	if(!is.environment(envir)) envir <- environment(as.formula(formula(x)))
	
	if(RE.keep) {
		if(inherits(x, c("mer", "merMod", "coxme", "lmekin"))) {
			cl$formula <- .nullREForm(as.formula(cl$formula))
			environment(cl$formula) <- envir
		} else if(inherits(x, "gamm")) {
			mefm <- x[[if("lme" %in% names(x)) "lme" else "mer"]]
			
			if(inherits(mefm, "merMod")) {
				Fun <- if(inherits(mefm, "glmerMod"))
					"glmer" else if(inherits(mefm, "lmerMod")) {
					cl$family <- NULL
					"lmer"
				}
				cl$REML <- as.logical(x$mer@devcomp$dims[['REML']])
				frm <- cl$formula
				frm[[3L]] <- call("+", 1, as.formula(cl$random)[[2L]])
				cl$random <- NULL
				environment(cl$formula) <- envir
			} else if (inherits(mefm, "lme")) {
				Fun <- "lme"
				cl$fixed <- update.formula(as.formula(cl$formula), . ~ 1)
				cl$formula <- cl$family <- NULL
				cl$method <- x$lme$method
				environment(cl$fixed) <- envir
			}
			cl[[1L]] <- as.symbol(Fun)
		} else if(inherits(x, c("glmmML", "glimML"))) {
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
			environment(cl$formula) <- envir
		} else if(inherits(x, "lme")) {
			cl$fixed <- update.formula(as.formula(cl$fixed), . ~ 1)
			environment(cl$fixed) <- envir
		} else {
			stop("do not know (yet) how to construct a null model with RE for class ",
				 prettyEnumStr(class(x), sep.last = ", "))					
		}
		return(if(evaluate) eval(cl, envir = envir) else cl)
	}	
	
	mClasses <- c("glmmML", "lm", "lme", "gls", "mer", "merMod", "lmekin",
				  "unmarkedFit", "coxph", "coxme", "zeroinfl", "gamm")
	mClass <- mClasses[inherits(x, mClasses, which = TRUE) != 0L][1]

	if(is.na(mClass)) mClass <- "default"
	formulaArgName <- "formula"
	Fun <- "glm"
	call2arg <- function(x) formals(match.fun(x[[1L]]))
	switch(mClass,
		glmmML = {
			if(is.null(cl$family)) cl$family <- as.name("binomial")
		}, gls = {
			formulaArgName <- "model"
			cl$weights <- NULL
		}, lme = {
			formulaArgName <- "fixed"
			cl$weights <- NULL
		}, lmekin =, merMod =, mer = {
			arg <- formals(match.fun(cl[[1L]]))
		}, unmarkedFit = {
			nm <- names(cl)[-1L]
			if("formula" %in% nm) {
				cl$formula <- ~1~1
			} else {
				formula.arg <- nm[grep(".+formula$", nm[1:7])]
				for (i in formula.arg) cl[[i]] <- ~1
			}
			cl$starts <- NULL
			Fun <- NA
		}, coxph =, coxme = {
			Fun <- "coxph"
			cl$formula <- update.formula(eval(cl$formula), . ~ 1)
		}, zeroinfl =, lm = {
			Fun <- NA
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
		}, gamm = {
			Fun <- "gam"
			cl$formula <- update.formula(as.formula(cl$formula), . ~ 1)
			cl$random <- NULL
		}, {
			stop("do not know (yet) how to construct a null model for class ",
				sQuote(class(x)))
		}
	)
	

	if(!is.na(Fun)) cl[[1L]] <- as.name(Fun)
	if(identical(Fun, "glm")) {
		if(formulaArgName != "formula")
			names(cl)[names(cl) == formulaArgName] <- "formula"
		cl$formula <- update(as.formula(cl$formula), . ~ 1)
		cl$method <- cl$start <- cl$offset <- contrasts <- NULL
	}
	cl <- cl[c(TRUE, names(cl)[-1L] %in% names(call2arg(cl)))]
	if(evaluate) eval(cl, envir = envir) else cl
}
		

# from lme4:::findbars:
.findbars <- function (term) {
    if (is.name(term) || !is.language(term))
        return(NULL)
    if (term[[1L]] == as.name("("))
        return(.findbars(term[[2L]]))
    if (!is.call(term))
        stop("term must be of class call")
    if (term[[1L]] == as.name("|"))
        return(term)
    if (length(term) == 2L)
        return(.findbars(term[[2L]]))
    c(.findbars(term[[2L]]), .findbars(term[[3L]]))
}

`.nullREForm` <-
function(formula) {
	re <- lapply(.findbars(formula), function(x) call("(", x))
	f <- 1
	for(i in seq_along(re)) f <- call("+", f, re[[i]])
	formula[[length(formula)]] <- f
	formula
}

`r.squaredLR` <-
function(x, null = NULL, null.RE = FALSE) {
	if(!missing(null) && !missing(null.RE))
		warning("argument 'null.RE' disregarded if 'null' is provided")
	if(is.null(null))
		null <- null.fit(x, TRUE, null.RE, parent.frame())

	#print(get_call(null))
	L0 <- as.vector(if(inherits(null, "glm")) logLik(null) else logLik(null, REML = FALSE))
	L1 <- if(inherits(x, "glm")) logLik(x) else logLik(x, REML = FALSE)
	n <- if(is.null(attr(L1, "nobs"))) nobs(x) else attr(L1, "nobs")
	ret <- 1 - exp(-2 / n * (as.vector(L1) - L0))
	max.r2 <- 1 - exp(2 / n * L0)
	attr(ret, "adj.r.squared") <- ret / max.r2
	ret
}

