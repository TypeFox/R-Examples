#-------------------------------------------------------------------------------
# Revision history:
# 2009-01-05: bug fix in Anova.II.lm(). J. Fox
# 2009-01-16: Cox models with clusters now handled. J. Fox
# 2009-09-16: reworked glm and lm methods to handle aliased parameters. J. Fox
# 2009-09-30: renamed "Anova" to "Analysis of Deviance" in output for some methods. J. Fox
# 2009-12-22: modified Anova.mlm() to handle a user-supplied within-subject model matrix. J. Fox
# 2009-12-28: named the components of P in Anova.III.mlm(). John
# 2010-01-01: Anova.II.mlm() now hands off (again) to Anova.III.mlm() when there
#             is only an intercept in the between-subjects model
# 2010-02-17: Fixed bug that caused some models with aliased coefficients to fail. J. Fox
# 2010-06-14: added wcrossprod and allow use of observation weights in Anova.mlm()
# 2010-06-28: Fixed Anova() tables for coxph and survreg models 
#             (failed because of changes in survival package.
# 2011-01-21: Added functions for mixed models. J. Fox
# 2011-01-25: Fixed Anova.polr() and Anova.multinom() to work with models with only one term. J. Fox
# 2011-05-19: local fixef() to avoid nlme/lme4 issues. J. Fox
# 2011-05-11: changed order of columns in ANOVA tables for mixed models. J. Fox
# 2011-11-27: added Anova.svyglm(). J. Fox
# 2011-12-31: fixed bug in Anova.II(and III).F.glm() when na.exclude used. J. Fox
# 2012-02-28: added test.statistic argument to Anova.mer(). J.Fox
# 2012-03-02: fixed test abbreviation of test.statistic argument to Anova.default()
#             called by other Anova() methods. J. Fox
# 2013-06-17: modified summary.Anova.mlm(), introduced print.summary.Anova.mlm(),
#             adapting code contributed by Gabriel Baud-Bovy. J. Fox
# 2013-06-20: added Anova.merMod() method. J. Fox
# 2013-06-22: tweaks to local fixef(). J. Fox
# 2013-06-22: test argument uniformly uses "Chisq" rather than "chisq". J. Fox
# 2013-08-19: replaced calls to print.anova(). J. Fox
# 2014-08-17: added calls to requireNamespace() and :: where needed (doesn't work for pbkrtest). J. Fox
# 2014-08-18: fixed bugs in Anova.survreg() for types II, III LR tests and Wald tests. J. Fox
# 2014-09-23: added Anova.rlm(). J. Fox
# 2014-10-10: removed MASS:: from calls to polr(). John
# 2014-12-18: check that residual df and SS are nonzero in Anova.lm(). John
# 2015-01-27: vcovAdj() and methods now imported from pbkrtest. John
# 2015-02-18: force evaluation of vcov. when it's a function. John
# 2015-04-30: don't allow error.estimate="dispersion" for F-tests in binomial
#             and Poission GLMs. John
# 2015-08-29: fixed Anova() for coxph models with clusters. John
# 2015-09-04: added support for coxme models. John
# 2015-09-11: modified Anova.default() to work with vglm objects from VGAM. John
# 2015-09-15: fixed Anova.default() so that F-tests work again. John
# 2015-11-13: modify Anova.coxph() to take account of method/ties argument. John
#-------------------------------------------------------------------------------

# Type II and III tests for linear, generalized linear, and other models (J. Fox)

ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
	# This function by Georges Monette
	# finds the conjugate complement of the proj of X in span(Z) wrt
	#    inner product ip
	# - assumes Z is of full column rank
	# - projects X conjugately wrt ip into span Z
	xq <- qr(t(Z) %*% ip %*% X)
	if (xq$rank == 0) return(Z)
	Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
}

relatives <- function(term, names, factors){
	is.relative <- function(term1, term2) {
		all(!(factors[,term1]&(!factors[,term2])))
	}
	if(length(names) == 1) return(NULL)
	which.term <- which(term==names)
	(1:length(names))[-which.term][sapply(names[-which.term], 
					function(term2) is.relative(term, term2))]
}


Anova <- function(mod, ...){
	UseMethod("Anova", mod)
}

# linear models

Anova.lm <- function(mod, error, type=c("II","III", 2, 3), 
		white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
        vcov.=NULL, singular.ok, ...){
    if (is.function(vcov.)) vcov. <- vcov.(mod)
    if (df.residual(mod) == 0) stop("residual df = 0")
    if (deviance(mod) < sqrt(.Machine$double.eps)) stop("residual sum of squares is 0 (within rounding error)")
	type <- as.character(type)
	white.adjust <- as.character(white.adjust)
	type <- match.arg(type)
	white.adjust <- match.arg(white.adjust)
	if (missing(singular.ok)){
		singular.ok <- type == "2" || type == "II"
	}
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	if (white.adjust != "FALSE"){
		if (white.adjust == "TRUE") white.adjust <- "hc3" 
		return(Anova.default(mod, type=type, vcov.=hccm(mod, type=white.adjust), test.statistic="F", 
						singular.ok=singular.ok, ...))
	}
    else if (!is.null(vcov.)) return(Anova.default(mod, type=type, vcov.=vcov., test.statistic="F", 
        singular.ok=singular.ok, ...))
	switch(type,
			II=Anova.II.lm(mod, error, singular.ok=singular.ok, ...),
			III=Anova.III.lm(mod, error, singular.ok=singular.ok, ...),
			"2"=Anova.II.lm(mod, error, singular.ok=singular.ok, ...),
			"3"=Anova.III.lm(mod, error, singular.ok=singular.ok,...))
}

Anova.aov <- function(mod, ...){
	class(mod) <- "lm"
	Anova.lm(mod, ...)
}

Anova.II.lm <- function(mod, error, singular.ok=TRUE, ...){
	if (!missing(error)){
		sumry <- summary(error, corr=FALSE)
		s2 <- sumry$sigma^2
		error.df <- error$df.residual
		error.SS <- s2*error.df
	}
	SS.term <- function(term){
		which.term <- which(term == names)
		subs.term <- which(assign == which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign == relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod)))
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(SS=NA, df=0))
		lh <- linearHypothesis(mod, hyp.matrix.term, 
				singular.ok=singular.ok, ...)
		abs(c(SS=lh$"Sum of Sq"[2], df=lh$Df[2]))
	}
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	assign <- mod$assign
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <-names[-1]
	n.terms <- length(names)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	sumry <- summary(mod, corr = FALSE)
	SS[n.terms + 1] <- if (missing(error)) sumry$sigma^2*mod$df.residual 
			else error.SS   
	df[n.terms + 1] <- if (missing(error)) mod$df.residual else error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	for (i in 1:n.terms){
		ss <- SS.term(names[i])
		SS[i] <- ss["SS"]
		df[i] <- ss["df"]
		f[i] <- df[n.terms+1]*SS[i]/(df[i]*SS[n.terms + 1])
		p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
	}    
	result <- data.frame(SS, df, f, p)
	row.names(result) <- c(names,"Residuals")
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

# type III

Anova.III.lm <- function(mod, error, singular.ok=FALSE, ...){
	if (!missing(error)){
		error.df <- df.residual(error)
		error.SS <- deviance(error)
	}
	else {
		error.df <- df.residual(mod)
		error.SS <- deviance(mod)
	}
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	Source <- term.names(mod)
	n.terms <- length(Source)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	assign <- mod$assign
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix) == 0){
			SS[term] <- NA
			df[term] <- 0
			f[term] <- NA
			p[term] <- NA
		}
		else {
			test <- if (missing(error)) linearHypothesis(mod, hyp.matrix, 
								singular.ok=singular.ok, ...)
					else linearHypothesis(mod, hyp.matrix, error.SS=error.SS, error.df=error.df, 
								singular.ok=singular.ok, ...)
			SS[term] <- test$"Sum of Sq"[2]
			df[term] <- test$"Df"[2]
			f[term] <- test$"F"[2]
			p[term] <- test$"Pr(>F)"[2]
		}
	}
	Source[n.terms + 1] <- "Residuals"
	SS[n.terms + 1] <- error.SS
	df[n.terms + 1] <- error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	result <- data.frame(SS, df, f, p)
	row.names(result) <- Source
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}

# generalized linear models

Anova.glm <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald", "F"), 
		error, error.estimate=c("pearson", "dispersion", "deviance"), singular.ok, ...){
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	if (missing(singular.ok)){
		singular.ok <- type == "2" || type == "II"
	}
	test.statistic <- match.arg(test.statistic)
	error.estimate <- match.arg(error.estimate)
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="II", singular.ok=singular.ok),
					F=Anova.II.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			III=switch(test.statistic,
					LR=Anova.III.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="III", singular.ok=singular.ok),
					F=Anova.III.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="II", singular.ok=singular.ok),
					F=Anova.II.F.glm(mod, error, error.estimate, singular.ok=singular.ok)),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.glm(mod, singular.ok=singular.ok),
					Wald=Anova.default(mod, type="III", singular.ok=singular.ok),
					F=Anova.III.F.glm(mod, error, error.estimate, singular.ok=singular.ok)))
}


# type III

# LR test

Anova.III.LR.glm <- function(mod, singular.ok=FALSE, ...){
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	Source <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(Source)
	p <- df <- LR <- rep(0, n.terms)
	dispersion <- summary(mod, corr = FALSE)$dispersion
	deviance <- deviance(mod)/dispersion
	for (term in 1:n.terms){
		mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
		LR[term] <- (mod.1$Deviance[2]/dispersion)-deviance
		df[term] <- mod.1$Df[2]
		p[term] <- pchisq(LR[term], df[term], lower.tail = FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- Source
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova","data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}

# F test

Anova.III.F.glm <- function(mod, error, error.estimate, singular.ok=FALSE, ...){
    if (!singular.ok && any(is.na(coef(mod))))
        stop("there are aliased coefficients in the model")
    fam <- family(mod)$family
    if ((fam == "binomial" || fam == "poisson") && error.estimate == "dispersion"){
        warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
        error.estimate <- "pearson"
    }
    if (missing(error)) error <- mod
    df.res <- df.residual(error)
    error.SS <- switch(error.estimate,
        pearson=sum(residuals(error, "pearson")^2, na.rm=TRUE),
        dispersion=df.res*summary(error, corr = FALSE)$dispersion,
        deviance=deviance(error))
    Source <- if (has.intercept(mod)) term.names(mod)[-1]
    else term.names(mod)
    n.terms <- length(Source)
    p <- df <- f <- SS <-rep(0, n.terms+1)
    f[n.terms+1] <- p[n.terms+1] <- NA
    df[n.terms+1] <- df.res
    SS[n.terms+1] <- error.SS
    dispersion <- error.SS/df.res
    deviance <- deviance(mod)
    for (term in 1:n.terms){
        mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
        df[term] <- mod.1$Df[2]
        SS[term] <- mod.1$Deviance[2] - deviance
        f[term] <- (SS[term]/df[term])/dispersion
        p[term] <- pf(f[term], df[term], df.res, lower.tail = FALSE)
    }
    result <- data.frame(SS, df, f, p)
    row.names(result) <- c(Source, "Residuals")
    names(result) <- c("SS", "Df", "F", "Pr(>F)")
    class(result) <- c("anova","data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
        paste("Response:", responseName(mod)), 
        paste("Error estimate based on",
            switch(error.estimate, 
                pearson="Pearson residuals", dispersion="estimated dispersion", 
                deviance="deviance"), "\n"))
    result
}

# type II

# LR test

Anova.II.LR.glm <- function(mod, singular.ok=TRUE, ...){
	if (!singular.ok && any(is.na(coef(mod))))
		stop("there are aliased coefficients in the model")
	# (some code adapted from drop1.glm)
	which.nms <- function(name) which(asgn == which(names == name))
	fac <- attr(terms(mod), "factors")
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- mod$y
	if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
	wt <- mod$prior.weights
	if (is.null(wt)) wt <- rep(1, length(y))
	asgn <- attr(X, 'assign')
	df <- p <- LR <- rep(0, n.terms)
	dispersion <- summary(mod, corr = FALSE)$dispersion
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
				family = mod$family, control = mod$control)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
							family = mod$family, control = mod$control)
				}
		dev.2 <- deviance(mod.2)
		df[term] <- df.residual(mod.1) - df.residual(mod.2)
		if (df[term] == 0) LR[term] <- p[term] <- NA
		else {
			LR[term] <- (dev.1 - dev.2)/dispersion
			p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
		}
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- 
			c("Analysis of Deviance Table (Type II tests)\n", paste("Response:", responseName(mod)))
	result
}


# F test

Anova.II.F.glm <- function(mod, error, error.estimate, singular.ok=TRUE, ...){
    # (some code adapted from drop1.glm)
    if (!singular.ok && any(is.na(coef(mod))))
        stop("there are aliased coefficients in the model")
    fam <- family(mod)$family
    if ((fam == "binomial" || fam == "poisson") && error.estimate == "dispersion"){
        warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
        error.estimate <- "pearson"
    }
    which.nms <- function(name) which(asgn == which(names == name))
    if (missing(error)) error <- mod
    df.res <- df.residual(error)
    error.SS <- switch(error.estimate,
        pearson = sum(residuals(error, "pearson")^2, na.rm=TRUE),
        dispersion = df.res*summary(error, corr = FALSE)$dispersion,
        deviance = deviance(error))
    fac <- attr(terms(mod), "factors")
    names <- if (has.intercept(mod)) term.names(mod)[-1]
    else term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- mod$y
    if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
    wt <- mod$prior.weights
    if (is.null(wt)) wt <- rep(1, length(y))
    asgn <- attr(X, 'assign')
    p <- df <- f <- SS <- rep(0, n.terms+1)
    f[n.terms+1] <- p[n.terms+1] <- NA
    df[n.terms+1] <- df.res
    SS[n.terms+1] <- error.SS
    dispersion <- error.SS/df.res
    for (term in 1:n.terms){
        rels <- names[relatives(names[term], names, fac)]
        exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
        mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
            family = mod$family, control = mod$control)
        dev.1 <- deviance(mod.1)
        mod.2 <- if (length(rels) == 0) mod
        else {
            exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
            glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
                family = mod$family, control = mod$control)
        }
        dev.2 <- deviance(mod.2)
        df[term] <- df.residual(mod.1) - df.residual(mod.2)
        if (df[term] == 0) SS[term] <- f[term] <- p[term] <- NA
        else {
            SS[term] <- dev.1 - dev.2
            f[term] <- SS[term]/(dispersion*df[term])
            p[term] <- pf(f[term], df[term], df.res, lower.tail=FALSE)
        }
    }
    result <- data.frame(SS, df, f, p)
    row.names(result) <- c(names, "Residuals")
    names(result) <- c("SS", "Df", "F", "Pr(>F)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
        paste("Response:", responseName(mod)),
        paste("Error estimate based on",
            switch(error.estimate,
                pearson="Pearson residuals", 
                dispersion="estimated dispersion", 
                deviance="deviance"), "\n"))
    result
}

# multinomial logit models (via multinom in the nnet package)

Anova.multinom <-
		function (mod, type = c("II", "III", 2, 3), ...)
{
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	switch(type,
			II = Anova.II.multinom(mod, ...),
			III = Anova.III.multinom(mod, ...),
			"2" = Anova.II.multinom(mod, ...),
			"3" = Anova.III.multinom(mod, ...))
}

Anova.II.multinom <- function (mod, ...)
{
	which.nms <- function(name) which(asgn == which(names ==
								name))
	fac <- attr(terms(mod), "factors")
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- mod$weights
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms) {
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
								which.nms)))
		mod.1 <-if (n.terms > 1) multinom(y ~ X[, -c(1, exclude.1)], weights=wt, trace=FALSE)
				else multinom(y ~ 1, weights=wt, race=FALSE)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0)
					mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					multinom(y ~ X[, -c(1, exclude.2)], weights=wt, trace=FALSE)
				}
		dev.2 <- deviance(mod.2)
		LR[term] <- dev.1 - dev.2
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
			paste("Response:", responseName(mod)))
	result
}

Anova.III.multinom <- function (mod, ...)
{
	names <- if (has.intercept(mod)) term.names(mod)[-1]
			else term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- mod$weights
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	deviance <- deviance(mod)
	for (term in 1:n.terms) {
		mod.1 <- if (n.terms > 1) multinom(y ~ X[, term != asgn][, -1], weights=wt, trace=FALSE)
				else multinom(y ~ 1, weights=wt, trace=FALSE)
		LR[term] <- deviance(mod.1) - deviance
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
			paste("Response:", responseName(mod)))
	result
}


# proportional-odds logit models (via polr in the MASS package)

Anova.polr <- function (mod, type = c("II", "III", 2, 3), ...)
{
	type <- as.character(type)
	type <- match.arg(type)
	if (has.intercept(mod) && length(coef(mod)) == 1 
			&& (type == "2" || type == "II")) {
		type <- "III"
		warning("the model contains only an intercept: Type III test substituted")
	}
	switch(type,
			II = Anova.II.polr(mod, ...),
			III = Anova.III.polr(mod, ...),
			"2" = Anova.II.polr(mod, ...),
			"3" = Anova.III.polr(mod, ...))
}

Anova.II.polr <- function (mod, ...)
{
  if (!requireNamespace("MASS")) stop("MASS package is missing")
	which.nms <- function(name) which(asgn == which(names ==
								name))
	fac <- attr(terms(mod), "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- model.weights(model.frame(mod))
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms) {
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
								which.nms)))
		mod.1 <- if (n.terms > 1) polr(y ~ X[, -c(1, exclude.1)], weights=wt)
				else polr(y ~ 1, weights=wt)
		dev.1 <- deviance(mod.1)
		mod.2 <- if (length(rels) == 0)
					mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					polr(y ~ X[, -c(1, exclude.2)], weights=wt)
				}
		dev.2 <- deviance(mod.2)
		LR[term] <- dev.1 - dev.2
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
			paste("Response:", responseName(mod)))
	result
}

Anova.III.polr <- function (mod, ...)
{
  if (!requireNamespace("MASS")) stop("MASS package is missing")
	names <- term.names(mod)
	n.terms <- length(names)
	X <- model.matrix(mod)
	y <- model.response(model.frame(mod))
	wt <- model.weights(model.frame(mod))
	asgn <- attr(X, "assign")
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	deviance <- deviance(mod)
	for (term in 1:n.terms) {
		mod.1 <- if (n.terms > 1) polr(y ~ X[, term != asgn][, -1], weights=wt)
				else polr(y ~ 1, weights=wt)
		LR[term] <- deviance(mod.1) - deviance
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
			paste("Response:", responseName(mod)))
	result
}

# multivariate linear models

# the following 3 functions copied from the stats package (not exported from stats)

Pillai <- function (eig, q, df.res) {
    test <- sum(eig/(1 + eig))
    p <- length(eig)
    s <- min(p, q)
    n <- 0.5 * (df.res - p - 1)
    m <- 0.5 * (abs(p - q) - 1)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * n + s + 1
    c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res) {
    test <- prod(1/(1 + eig))
    p <- length(eig)
    tmp1 <- df.res - 0.5 * (p - q + 1)
    tmp2 <- (p * q - 2)/4
    tmp3 <- p^2 + q^2 - 5
    tmp3 <- if (tmp3 > 0) 
        sqrt(((p * q)^2 - 4)/tmp3)
    else 1
    c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
      p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res) {
    test <- sum(eig)
    p <- length(eig)
    m <- 0.5 * (abs(p - q) - 1)
    n <- 0.5 * (df.res - p - 1)
    s <- min(p, q)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * (s * n + 1)
    c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res) {
    p <- length(eig)
    test <- max(eig)
    tmp1 <- max(p, q)
    tmp2 <- df.res - tmp1 + q
    c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}

has.intercept.mlm <- function (model, ...) 
	any(row.names(coefficients(model)) == "(Intercept)")


Anova.mlm <- function(mod, type=c("II","III", 2, 3), SSPE, error.df, idata, 
		idesign, icontrasts=c("contr.sum", "contr.poly"), imatrix,
		test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...){
	wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (missing(SSPE)) SSPE <- wcrossprod(residuals(mod), w=wts)
	if (missing(idata)) {
		idata <- NULL
		idesign <- NULL
	}
	if (missing(imatrix)) imatrix <- NULL
	error.df <- if (missing(error.df)) df.residual(mod)
			else error.df
	switch(type,
			II=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			III=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			"2"=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
			"3"=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...))
}

Anova.III.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
	intercept <- has.intercept(mod)
	V <- solve(crossprod(model.matrix(mod)))
	p <- nrow(coefficients(mod))
	I.p <- diag(p)
	terms <- term.names(mod)
	n.terms <- length(terms)
	assign <- mod$assign
	if (is.null(idata) && is.null(imatrix)){
		if ((n.terms == 0) && intercept) {
			Test <- linearHypothesis(mod, 1, SSPE=SSPE, ...)
			result <- list(SSP=Test$SSPH, SSPE=SSPE, df=1, error.df=error.df,
					terms="(Intercept)", repeated=FALSE, type="III", test=test)
			class(result) <- "Anova.mlm"
			return(result)
		}
		SSP <- as.list(rep(0, n.terms))
		df <- rep(0, n.terms)
		names(df) <- names(SSP) <- terms
		for (term in 1:n.terms){
			subs <- which(assign == term - intercept)
			hyp.matrix <- I.p[subs,,drop=FALSE]
			Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, ...)
			SSP[[term]] <- Test$SSPH
			df[term]<- length(subs)
		}
		result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
				repeated=FALSE, type="III", test=test)
	}
	else {
		if (!is.null(imatrix)){
			X.design <- do.call(cbind, imatrix)
			ncols <- sapply(imatrix, ncol)
			end <- cumsum(ncols)
			start <- c(1, (end + 1))[-(length(end) + 1)]
			cols <- mapply(seq, from=start, to=end)
			iterms <- names(end)
			names(cols) <- iterms
			check.imatrix(X.design, iterms)
		}
		else {
			if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
			for (i in 1:length(idata)){
				if (is.null(attr(idata[,i], "contrasts"))){
					contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
							else icontrasts[1]
				}
			}
			X.design <- model.matrix(idesign, data=idata)
			i.intercept <- has.intercept(X.design)
			iterms <- term.names(idesign)
			if (i.intercept) iterms <- c("(Intercept)", iterms)
			check.imatrix(X.design)
		}
		df <- rep(0, n.terms*length(iterms))
		hnames <- rep("", length(df))
		P <- SSPEH <- SSP <- as.list(df)
		i <- 0
		for (iterm in iterms){
			for (term in 1:n.terms){
				subs <- which(assign == term - intercept)
				hyp.matrix <- I.p[subs,,drop=FALSE]
				i <- i + 1
				Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, 
						idata=idata, idesign=idesign, icontrasts=icontrasts, iterms=iterm, 
						check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)
				SSP[[i]] <- Test$SSPH
				SSPEH[[i]] <- Test$SSPE
				P[[i]] <- Test$P
				df[i] <- length(subs)
				hnames[i] <- if (iterm == "(Intercept)") terms[term]
						else if (terms[term] == "(Intercept)") iterm
						else paste(terms[term], ":", iterm, sep="")
			}
		}
		names(df) <- names(SSP) <- names(SSPEH) <- names(P) <- hnames
		result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
				terms=hnames, repeated=TRUE, type="III", test=test, 
				idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix,
				singular=Test$singular)       
	}
	class(result) <- "Anova.mlm"
	result
}

Anova.II.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
	wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
	V <- solve(wcrossprod(model.matrix(mod), w=wts))
	SSP.term <- function(term, iterm){
		which.term <- which(term == terms)
		subs.term <- which(assign == which.term)
		relatives <- relatives(term, terms, fac)
		subs.relatives <- NULL
		for (relative in relatives) subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),,drop=FALSE]
		if (missing(iterm)){
			SSP1 <- if (length(subs.relatives) == 0) 0 
					else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, singular.ok=TRUE, ...)$SSPH
			SSP2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, singular.ok=TRUE, ...)$SSPH
			return(SSP2 - SSP1)
		}
		else {
			SSP1 <- if (length(subs.relatives) == 0) 0 
					else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
								idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], singular.ok=TRUE, ...)$SSPH
			lh2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, 
					idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], singular.ok=TRUE, ...)
			return(list(SSP = lh2$SSPH - SSP1, SSPE=lh2$SSPE, P=lh2$P, singular=lh2$singular))
		}
	}
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- nrow(coefficients(mod))
	I.p <- diag(p)
	assign <- mod$assign
	terms <- term.names(mod)
	if (intercept) terms <- terms[-1]
	n.terms <- length(terms)
	if (n.terms == 0){
		message("Note: model has only an intercept; equivalent type-III tests substituted.")
		return(Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...))
	}
	if (is.null(idata) && is.null(imatrix)){
		SSP <- as.list(rep(0, n.terms))
		df <- rep(0, n.terms)
		names(df) <- names(SSP) <- terms
		for (i in 1:n.terms){
			SSP[[i]] <- SSP.term(terms[i])
			df[i]<- df.terms(mod, terms[i])
		}    
		result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
				repeated=FALSE, type="II", test=test)
	}
	else {
		if (!is.null(imatrix)){
			X.design <- do.call(cbind, imatrix)
			ncols <- sapply(imatrix, ncol)
			end <- cumsum(ncols)
			start <- c(1, (end + 1))[-(length(end) + 1)]
			cols <- mapply(seq, from=start, to=end)
			iterms <- names(end)
			names(cols) <- iterms
			check.imatrix(X.design, iterms)     
		}
		else {
			if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
			for (i in 1:length(idata)){
				if (is.null(attr(idata[,i], "contrasts"))){
					contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
							else icontrasts[1]
				}
			}
			X.design <- model.matrix(idesign, data=idata)
			iintercept <- has.intercept(X.design)
			iterms <- term.names(idesign)
			if (iintercept) iterms <- c("(Intercept)", iterms)
			check.imatrix(X.design)
		}
		df <- rep(0, (n.terms + intercept)*length(iterms))
		hnames <- rep("", length(df))
		P <- SSPEH <- SSP <- as.list(df)
		i <- 0
		for (iterm in iterms){
			if (intercept){
				i <- i + 1
				hyp.matrix.1 <- I.p[-1,,drop=FALSE]
				SSP1 <- linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
						idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
						check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)$SSPH
				lh2 <- linearHypothesis(mod, I.p, SSPE=SSPE, V=V, 
						idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
						check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)
				SSP[[i]] <- lh2$SSPH - SSP1
				SSPEH[[i]] <- lh2$SSPE
				P[[i]] <- lh2$P
				df[i] <- 1
				hnames[i] <- iterm
			}
			for (term in 1:n.terms){
				subs <- which(assign == term)
				i <- i + 1
				Test <- SSP.term(terms[term], iterm)
				SSP[[i]] <- Test$SSP
				SSPEH[[i]] <- Test$SSPE
				P[[i]] <- Test$P
				df[i]<- length(subs)
				hnames[i] <- if (iterm == "(Intercept)") terms[term]
						else paste(terms[term], ":", iterm, sep="")
			}
		}
		names(df) <- names(P) <- names(SSP) <- names(SSPEH) <- hnames
		result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
				terms=hnames, repeated=TRUE, type="II", test=test,
				idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix,
				singular=Test$singular)       
	}
	class(result) <- "Anova.mlm"
	result
}

print.Anova.mlm <- function(x, ...){
	if ((!is.null(x$singular)) && x$singular) stop("singular error SSP matrix; multivariate tests unavailable\ntry summary(object, multivariate=FALSE)")
	test <- x$test
	repeated <- x$repeated
	ntests <- length(x$terms)
	tests <- matrix(NA, ntests, 4)
	if (!repeated) SSPE.qr <- qr(x$SSPE) 
	for (term in 1:ntests){
		# some of the code here adapted from stats:::summary.manova
		eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr,
								x$SSP[[term]]), symmetric = FALSE)$values)
		tests[term, 1:4] <- switch(test,
				Pillai = Pillai(eigs, x$df[term], x$error.df),
				Wilks = Wilks(eigs, x$df[term], x$error.df),
				"Hotelling-Lawley" = HL(eigs, x$df[term], x$error.df),
				Roy = Roy(eigs, x$df[term], x$error.df))
	}
	ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
	ok <- !is.na(ok) & ok
	tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
					lower.tail = FALSE))
	rownames(tests) <- x$terms
	colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
	tests <- structure(as.data.frame(tests), 
			heading = paste("\nType ", x$type, if (repeated) " Repeated Measures",
					" MANOVA Tests: ", test, " test statistic", sep=""), 
			class = c("anova", "data.frame"))
	print(tests)      
	invisible(x)
}

# path <-  "D:/R-package-sources/car/R"
# files <- list.files(path, pattern=".*\\.R")
# files <- paste(path, files, sep="/")
# for (file in files) source(file)

# summary.Anova.mlm and print.summary.Anova.mlm methods
#  with contributions from Gabriel Baud-Bovy
summary.Anova.mlm <- function (object, test.statistic, univariate=TRUE, multivariate=TRUE, ...) {
    GG <- function(SSPE, P) { # Greenhouse-Geisser correction
        p <- nrow(SSPE)
        if (p < 2) 
            return(NA)
        lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
        lambda <- lambda[lambda > 0]
        ((sum(lambda)/p)^2)/(sum(lambda^2)/p)
    }
    HF <- function(gg, error.df, p) { # Huynh-Feldt correction
        ((error.df + 1) * p * gg - 2)/(p * (error.df - p * gg))
    }
    mauchly <- function(SSD, P, df) {
        # most of this function borrowed from stats:::mauchly.test.SSD
        if (nrow(SSD) < 2) 
            return(c(NA, NA))
        Tr <- function(X) sum(diag(X))
        p <- nrow(P)
        I <- diag(p)
        Psi <- t(P) %*% I %*% P
        B <- SSD
        pp <- nrow(SSD)
        U <- solve(Psi, B)
        n <- df
        logW <- log(det(U)) - pp * log(Tr(U/pp))
        rho <- 1 - (2 * pp^2 + pp + 2)/(6 * pp * n)
        w2 <- (pp + 2) * (pp - 1) * (pp - 2) * (2 * pp^3 + 6 * 
                pp^2 + 3 * p + 2)/(288 * (n * pp * rho)^2)
        z <- -n * rho * logW
        f <- pp * (pp + 1)/2 - 1
        Pr1 <- pchisq(z, f, lower.tail = FALSE)
        Pr2 <- pchisq(z, f + 4, lower.tail = FALSE)
        pval <- Pr1 + w2 * (Pr2 - Pr1)
        c(statistic = c(W = exp(logW)), p.value = pval)
    }
    if (missing(test.statistic)) 
        test.statistic <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
    test.statistic <- match.arg(test.statistic, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), several.ok = TRUE)
    nterms <- length(object$terms)
    summary.object <- list(type=object$type, repeated=object$repeated, 
        multivariate.tests=NULL, univariate.tests=NULL, 
        pval.adjustments=NULL, sphericity.tests=NULL)
    if (multivariate){
        summary.object$multivariate.tests <- vector(nterms, mode="list")
        names(summary.object$multivariate.tests) <- object$terms 
        summary.object$SSPE <- object$SSPE
        for (term in 1:nterms) {
            hyp <- list(SSPH = object$SSP[[term]], 
                SSPE = if (object$repeated) object$SSPE[[term]] else object$SSPE, 
                P = if (object$repeated) object$P[[term]] else NULL, 
                test = test.statistic, df = object$df[term], 
                df.residual = object$error.df, title = object$terms[term])
            class(hyp) <- "linearHypothesis.mlm"
            summary.object$multivariate.tests[[term]] <- hyp
        }
    }
    if (object$repeated && univariate) {
        singular <- object$singular
        error.df <- object$error.df
        table <- matrix(0, nterms, 6)
        table2 <- matrix(0, nterms, 4)
        table3 <- matrix(0, nterms, 2)
        rownames(table3) <- rownames(table2) <- rownames(table) <- object$terms
        colnames(table) <- c("SS", "num Df", "Error SS", "den Df", "F", "Pr(>F)")
        colnames(table2) <- c("GG eps", "Pr(>F[GG])", "HF eps","Pr(>F[HF])")
        colnames(table3) <- c("Test statistic", "p-value")
        if (singular) 
            warning("Singular error SSP matrix:\nnon-sphericity test and corrections not available")
        for (term in 1:nterms) {
            SSP <- object$SSP[[term]]
            SSPE <- object$SSPE[[term]]
            P <- object$P[[term]]
            p <- ncol(P)
            PtPinv <- solve(t(P) %*% P)
            gg <- if (!singular) GG(SSPE, P) else NA
            table[term, "SS"] <- sum(diag(SSP %*% PtPinv))
            table[term, "Error SS"] <- sum(diag(SSPE %*% PtPinv))
            table[term, "num Df"] <- object$df[term] * p
            table[term, "den Df"] <- error.df * p
            table[term, "F"] <- (table[term, "SS"]/table[term, "num Df"])/
                (table[term, "Error SS"]/table[term, "den Df"])
            table[term, "Pr(>F)"] <- pf(table[term, "F"], table[term, "num Df"], table[term, "den Df"], 
                lower.tail = FALSE)
            table2[term, "GG eps"] <- gg
            table2[term, "HF eps"] <- if (!singular) HF(gg, error.df, p) else NA
            table3[term, ] <- if (!singular) mauchly(SSPE, P, object$error.df) else NA
        }
        table3 <- na.omit(table3)
        if (nrow(table3) > 0) {
            table2[, "Pr(>F[GG])"] <- pf(table[, "F"], table2[, "GG eps"] * 
                    table[, "num Df"], table2[, "GG eps"] * table[, "den Df"], 
                lower.tail = FALSE)
            table2[, "Pr(>F[HF])"] <- pf(table[, "F"], pmin(1, table2[, "HF eps"]) * 
                    table[, "num Df"], pmin(1, table2[, "HF eps"]) * table[, "den Df"], 
                lower.tail = FALSE)
            table2 <- na.omit(table2)
            if (any(table2[, "HF eps"] > 1)) warning("HF eps > 1 treated as 1")
        }
        class(table3) <- class(table) <- "anova"
        summary.object$univariate.tests <- table
        summary.object$pval.adjustments <- table2
        summary.object$sphericity.tests <- table3
    }
    class(summary.object) <- "summary.Anova.mlm"
    summary.object
}

print.summary.Anova.mlm <- function(x, digits = getOption("digits"), ... ) {
    if (!is.null(x$multivariate.tests)) {
        cat(paste("\nType ", x$type, if (x$repeated) 
            " Repeated Measures", " MANOVA Tests:\n", sep = ""))
        if (!x$repeated) {
            cat("\nSum of squares and products for error:\n")
            print(x$SSPE, digits = digits)
        }
        for (term in 1:length(x$multivariate.tests)) {
            cat(paste("\n------------------------------------------\n", 
                "\nTerm:", names(x$multivariate.tests)[term], "\n"))
            print(x$multivariate.tests[[term]], digits = digits, SSPE = x$repeated, ...)
        }
    }
    if  (!is.null(x$univariate.tests)) {
        cat("\nUnivariate Type", x$type, "Repeated-Measures ANOVA Assuming Sphericity\n\n")
        print(x$univariate.tests)
        if (nrow(x$sphericity.tests) > 0) {
            cat("\n\nMauchly Tests for Sphericity\n\n")
            print(x$sphericity.tests)
            cat("\n\nGreenhouse-Geisser and Huynh-Feldt Corrections\n", 
                "for Departure from Sphericity\n\n")
            table <- x$pval.adjustments[, 1:2, drop = FALSE]
            class(table) <- "anova"
            print(table)
            cat("\n")
            table <- x$pval.adjustments[, 3:4, drop = FALSE]
            class(table)
            print(table)
        }
    }
    invisible(x)
}

Anova.manova <- function(mod, ...){
	class(mod) <- c("mlm", "lm")
	Anova(mod, ...)
}

Manova <- function(mod, ...){
	UseMethod("Manova")
}

Manova.mlm <- function(mod, ...){
	Anova(mod, ...)
}

# Cox regression models

df.residual.coxph <- function(object, ...){    
	object$n - sum(!is.na(coef(object)))
}

alias.coxph <- function(model){
	if(any(which <- is.na(coef(model)))) return(list(Complete=which))
	else list()
}

logLik.coxph <- function(object, ...) object$loglik[2]

Anova.coxph <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
		warning("LR tests unavailable with robust variances\nWald tests substituted")
		test.statistic <- "Wald"
	}
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.coxph(mod),
					Wald=Anova.default(mod, type="II", test.statistic="Chisq", vcov.=vcov(mod))),
			III=switch(test.statistic,
					LR=Anova.III.LR.coxph(mod),
					Wald=Anova.default(mod, type="III", test.statistic="Chisq", vcov.=vcov(mod))),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.coxph(mod),
					Wald=Anova.default(mod, type="II", test.statistic="Chisq", vcov.=vcov(mod))),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.coxph(mod),
					Wald=Anova.default(mod, type="III", test.statistic="Chisq", vcov.=vcov(mod))))
}

Anova.II.LR.coxph <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod, test="Chisq"))
	method <- mod$method
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	p <- LR <- rep(0, n.terms)
	df <- df.terms(mod)
	for (term in 1:n.terms){
		rels <- names[relatives(names[term], names, fac)]
		exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
		mod.1 <- survival::coxph(mod$y ~ X[, -exclude.1, drop = FALSE], method=method)
		loglik.1 <- logLik(mod.1)
		mod.2 <- if (length(rels) == 0) mod
				else {
					exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
					survival::coxph(mod$y ~ X[, -exclude.2, drop = FALSE], method=method)
				}
		loglik.2 <- logLik(mod.2)
		LR[term] <- -2*(loglik.1 - loglik.2)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
	result
}

Anova.III.LR.coxph <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
	which.nms <- function(name) which(asgn == which(names == name))
	fac <-attr(terms(mod), "factors")
	names <- term.names(mod)
	n.terms <- length(names)
	if (n.terms < 2) return(anova(mod, test="Chisq"))
	method <- mod$method
	X <- model.matrix(mod)
	asgn <- attr(X, 'assign')
	df <- df.terms(mod)
	LR <- p <- rep(0, n.terms)
	loglik1 <- logLik(mod)
	for (term in 1:n.terms){
		mod.0 <- survival::coxph(mod$y ~ X[, -which.nms(names[term])], method=method)
		LR[term] <- -2*(logLik(mod.0) - loglik1)
		p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
	}
	result <- data.frame(LR, df, p)
	row.names(result) <- names
	names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
	result
}

# parametric survival regression models

alias.survreg <- function(model){
	if(any(which <- diag(vcov(model)) < 1e-10)) return(list(Complete=which))
	else list()
}

logLik.survreg <- function(object, ...) object$loglik[2]

Anova.survreg <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
		warning("LR tests unavailable with robust variances\nWald tests substituted")
		test.statistic <- "Wald"
	}
	switch(type,
			II=switch(test.statistic,
					LR=Anova.II.LR.survreg(mod),
					Wald=Anova.II.Wald.survreg(mod)),
			III=switch(test.statistic,
					LR=Anova.III.LR.survreg(mod),
					Wald=Anova.III.Wald.survreg(mod)),
			"2"=switch(test.statistic,
					LR=Anova.II.LR.survreg(mod),
					Wald=Anova.II.Wald.survreg(mod)),
			"3"=switch(test.statistic,
					LR=Anova.III.LR.survreg(mod),
					Wald=Anova.III.Wald.survreg(mod)))
}

Anova.II.LR.survreg <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  dist <- mod$dist
  scale <- mod$call$scale
  weights <- model.frame(mod)$"(weights)"
  arg.list <- list(dist=dist)
  if (!is.null(scale)) arg.list$scale <- scale
  if (!is.null(weights)) arg.list$weights <- weights
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  asgn <- asgn[asgn != 0]
  if (has.intercept(mod)){
    int <- which(names == "(Intercept)")
    X <- X[, -int]
    names <- names[-int]
  }
  n.terms <- length(names)
  if (n.terms < 2) return(anova(mod))
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  y <- model.frame(mod)[,1]
  for (term in 1:n.terms){
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
    arg.list$formula <- y ~ X[, -exclude.1, drop = FALSE]
    mod.1 <- do.call(survival::survreg, arg.list)
    #    mod.1 <- survival::survreg(y ~ X[, -exclude.1, drop = FALSE])
    loglik.1 <- logLik(mod.1)
    mod.2 <- if (length(rels) == 0) mod
    else {
      arg.list$formula <- y ~ X[, -exclude.2, drop = FALSE]
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      do.call(survival::survreg, arg.list)     
      #     survival::survreg(y ~ X[, -exclude.2, drop = FALSE])
    }
    loglik.2 <- logLik(mod.2)
    LR[term] <- -2*(loglik.1 - loglik.2)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
  result
}

Anova.III.LR.survreg <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  dist <- mod$dist
  scale <- mod$call$scale
  weights <- model.frame(mod)$"(weights)"
  arg.list <- list(dist=dist)
  if (!is.null(scale)) arg.list$scale <- scale
  if (!is.null(weights)) arg.list$weights <- weights
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  asgn <- asgn[asgn != 0]
  if (has.intercept(mod)){
    int <- which(names == "(Intercept)")
    X <- X[, -int]
    names <- names[-int]
  }
  n.terms <- length(names)
  if (n.terms < 2) return(anova(mod))
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  y <- model.frame(mod)[,1]
  loglik1 <- logLik(mod)
  for (term in 1:n.terms){
    arg.list$formula <- y ~ X[, -which.nms(names[term])]
    mod.0 <- do.call(survival::survreg, arg.list)
    #    mod.0 <- survival::survreg(y ~ X[, -which.nms(names[term])])
    LR[term] <- -2*(logLik(mod.0) - loglik1)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
  result
}

Anova.II.Wald.survreg <- function(mod){
	V <- vcov(mod)
	p <- which(rownames(V) == "Log(scale)")
	if (length(p) > 0) V <- V[-p, -p]
	Anova.II.default(mod, V, test="Chisq")
}

Anova.III.Wald.survreg <- function(mod){
	V <- vcov(mod)
	p <- which(rownames(V) == "Log(scale)")
	if (length(p) > 0) V <- V[-p, -p]
	Anova.III.default(mod, V, test="Chisq")
}

# Default Anova() method: requires methods for vcov() (if vcov. argument not specified) and coef().

Anova.default <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Chisq", "F"), 
		vcov.=vcov(mod), singular.ok, ...){
    if (is.function(vcov.)) vcov. <- vcov.(mod)
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			III=Anova.III.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			"2"=Anova.II.default(mod, vcov., test.statistic, singular.ok=singular.ok),
			"3"=Anova.III.default(mod, vcov., test.statistic, singular.ok=singular.ok))
}

assignVector <- function(model){
    m <- model.matrix(model)
    assign <- attr(m, "assign")
    if (!is.null(assign)) return (assign)
    m <- model.matrix(formula(model), data=model.frame(model))
    assign <- attr(m, "assign")
    if (!has.intercept(model)) assign <- assign[assign != 0]
    assign
}

Anova.II.default <- function(mod, vcov., test, singular.ok=TRUE, ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- if (is.list(assign)) assign[[which.term]] else which(assign == which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives){
		    sr <- if (is.list(assign)) assign[[relative]] else which(assign == relative)
			subs.relatives <- c(subs.relatives, sr)
		}
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis.default(mod, hyp.matrix.term, 
				vcov.=vcov., test=test, singular.ok=singular.ok, ...)
		if (test=="Chisq") c(statistic=hyp$Chisq[2], df=hyp$Df[2])
		else c(statistic=hyp$F[2], df=hyp$Df[2])
	}
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	assign <- assignVector(mod) # attr(model.matrix(mod), "assign")
	if (!is.list(assign)) assign[!not.aliased] <- NA
	else if (intercept) assign <- assign[-1]
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	df <- c(rep(0, n.terms), df.residual(mod))
	if (inherits(mod, "coxph")){
		assign <- assign[assign != 0]
		clusters <- grep("^cluster\\(", names)
		if (length(clusters) > 0) {
			names <- names[-clusters]
			df <- df[-clusters]
			n.terms <- n.terms - length(clusters)
		}
	}
#	if (inherits(mod, "plm")) assign <- assign[assign != 0]
	p <- teststat <- rep(0, n.terms + 1)
	teststat[n.terms + 1] <- p[n.terms + 1] <- NA
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		p[i] <- if (test == "Chisq") 
					pchisq(teststat[i], df[i], lower.tail=FALSE) 
				else pf(teststat[i], df[i], df[n.terms + 1], lower.tail=FALSE)
	}
	result <- if (test == "Chisq"){ 
	    if (length(df) == n.terms + 1) df <- df[1:n.terms]
	    data.frame(df, teststat[!is.na(teststat)], p[!is.na(teststat)])
	}
	else data.frame(df, teststat, p)
	if (nrow(result) == length(names) + 1) names <- c(names,"Residuals")
	row.names(result) <- names
	names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
					else "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.III.default <- function(mod, vcov., test, singular.ok=FALSE, ...){
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- assignVector(mod) # attr(model.matrix(mod), "assign")
	df <- c(rep(0, n.terms), df.residual(mod))
	if (inherits(mod, "coxph")){
		if (intercept) names <- names[-1]
		assign <- assign[assign != 0]
		clusters <- grep("^cluster\\(", names)
		if (length(clusters) > 0) {
			names <- names[-clusters]
			df <- df[-clusters]
			n.terms <- n.terms - length(clusters)
		}
	}
#	if (inherits(mod, "plm")) assign <- assign[assign != 0]
	if (intercept) df[1] <- sum(grepl("^\\(Intercept\\)", names(coef(mod))))
	teststat <- rep(0, n.terms + 1)
	p <- rep(0, n.terms + 1)
	teststat[n.terms + 1] <- p[n.terms + 1] <- NA
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- if (is.list(assign)) assign[[term]] else which(assign == term - intercept)    
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis.default(mod, hyp.matrix, 
					vcov.=vcov., test=test, singular.ok=singular.ok, ...)
			teststat[term] <- if (test=="Chisq") hyp$Chisq[2] else hyp$F[2]
			df[term] <- abs(hyp$Df[2])
			p[term] <- if (test == "Chisq") 
						pchisq(teststat[term], df[term], lower.tail=FALSE) 
					else pf(teststat[term], df[term], df[n.terms + 1], lower.tail=FALSE)
		}
	}
	result <- if (test == "Chisq"){ 
	    if (length(df) == n.terms + 1) df <- df[1:n.terms]
	    data.frame(df, teststat[!is.na(teststat)], p[!is.na(teststat)])
	}
	else data.frame(df, teststat, p)
	if (nrow(result) == length(names) + 1) names <- c(names,"Residuals")
	row.names(result) <- names
	names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
					else "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

## functions for mixed models

# the following function, not exported, to make car consistent with CRAN and development versions of lme4 and with nlme

fixef <- function (object){
	if (isS4(object)) {
        if (!inherits(object, "merMod")) object@fixef
        else lme4::fixef(object)
	}
    else object$coefficients$fixed
}

Anova.merMod <- function(mod, type=c("II","III", 2, 3), 
                         test.statistic=c("Chisq", "F"),
                         vcov.=vcov(mod), singular.ok, ...){
    if (is.function(vcov.)) vcov. <- vcov.(mod)
    type <- as.character(type)
    type <- match.arg(type)
    test.statistic <- match.arg(test.statistic)
    if (missing(singular.ok))
        singular.ok <- type == "2" || type == "II"
    Anova.mer(mod=mod, type=type, test.statistic=test.statistic, vcov.=vcov.,
              singular.ok=singular.ok, ...)
}

Anova.mer <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Chisq", "F"),
		vcov.=vcov(mod), singular.ok, ...){
    if (is.function(vcov.)) vcov. <- vcov.(mod)
	type <- as.character(type)
	type <- match.arg(type)
	test.statistic <- match.arg(test.statistic)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.mer(mod, test=test.statistic, vcov., singular.ok=singular.ok),
			III=Anova.III.mer(mod, test=test.statistic,  vcov., singular.ok=singular.ok),
			"2"=Anova.II.mer(mod, test=test.statistic, vcov., singular.ok=singular.ok),
			"3"=Anova.III.mer(mod, test=test.statistic, vcov., singular.ok=singular.ok))
}

Anova.II.mer <- function(mod, vcov., singular.ok=TRUE, test=c("Chisq", "F"), ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- which(assign==which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis(mod, hyp.matrix.term, 
				vcov.=vcov., singular.ok=singular.ok, test=test, ...)
		if (test == "Chisq") return(c(statistic=hyp$Chisq[2], df=hyp$Df[2]))
		else return(c(statistic=hyp$F[2], df=hyp$Df[2], res.df=hyp$Res.Df[2]))
	}
	test <- match.arg(test)
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	if (!missing(vcov.)){
		vcov. <- if (test == "F"){
#					if (!require("pbkrtest")) stop("pbkrtest package required for F-tests on linear mixed model")
					as.matrix(vcovAdj(mod, details=0))
				}
				else vcov(mod)
	}
	assign <- attr(model.matrix(mod), "assign")
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	p <- teststat <- df <- res.df <- rep(0, n.terms)
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		res.df[i] <- hyp["res.df"]
		p[i] <- if (test == "Chisq") pchisq(teststat[i], df[i], lower.tail=FALSE) 
				else pf(teststat[i], df[i], res.df[i], lower.tail=FALSE)
	} 
	if (test=="Chisq"){
		result <- data.frame(teststat, df, p)
		row.names(result) <- names
		names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
		class(result) <- c("anova", "data.frame")
		attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald chisquare tests)\n", 
				paste("Response:", responseName(mod)))
	}
	else {
		result <- data.frame(teststat, df, res.df, p)
		row.names(result) <- names
		names(result) <- c ("F", "Df", "Df.res", "Pr(>F)")
		class(result) <- c("anova", "data.frame")
		attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)\n", 
				paste("Response:", responseName(mod)))
	}
	result
}

Anova.III.mer <- function(mod, vcov., singular.ok=FALSE, test=c("Chisq", "F"), ...){
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- attr(model.matrix(mod), "assign")
	p <- teststat <- df <- res.df <- rep(0, n.terms)
	if (intercept) df[1] <- 1
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	if (!missing(vcov.)){
		vcov. <- if (test == "F"){
#					if (!require("pbkrtest")) stop("pbkrtest package required for F-tests on linear mixed model")
					as.matrix(vcovAdj(mod, details=0))
				}
				else vcov(mod)
	}
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)        
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis(mod, hyp.matrix, test=test,
					vcov.=vcov., singular.ok=singular.ok, ...)
			if (test == "Chisq"){
				teststat[term] <-  hyp$Chisq[2] 
				df[term] <- abs(hyp$Df[2])
				p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
			}
			else{
				teststat[term] <-  hyp$F[2]
				df[term] <- abs(hyp$Df[2])
				res.df[term]=hyp$Res.Df[2]
				p[term] <- pf(teststat[term], df[term], res.df[term], lower.tail=FALSE)
			}
		}
	}
	if (test == "Chisq"){
		result <- data.frame(teststat, df, p)
		row.names(result) <- names
		names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
		class(result) <- c("anova", "data.frame")
		attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald chisquare tests)\n", 
				paste("Response:", responseName(mod)))
	}
	else {
		result <- data.frame(teststat, df, res.df, p)
		row.names(result) <- names
		names(result) <- c ("F", "Df", "Df.res", "Pr(>F)")
		class(result) <- c("anova", "data.frame")
		attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)\n", 
				paste("Response:", responseName(mod)))
	}
	result
}

Anova.lme <- function(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod), singular.ok, ...){
    if (is.function(vcov.)) vcov. <- vcov.(mod)
	type <- as.character(type)
	type <- match.arg(type)
	if (missing(singular.ok))
		singular.ok <- type == "2" || type == "II"
	switch(type,
			II=Anova.II.lme(mod, vcov., singular.ok=singular.ok),
			III=Anova.III.lme(mod, vcov., singular.ok=singular.ok),
			"2"=Anova.II.lme(mod, vcov., singular.ok=singular.ok),
			"3"=Anova.III.lme(mod, vcov., singular.ok=singular.ok))
}

Anova.II.lme <- function(mod, vcov., singular.ok=TRUE, ...){
	hyp.term <- function(term){
		which.term <- which(term==names)
		subs.term <- which(assign==which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign==relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(statistic=NA, df=0))            
		hyp <- linearHypothesis(mod, hyp.matrix.term, 
				vcov.=vcov., singular.ok=singular.ok, ...)
		c(statistic=hyp$Chisq[2], df=hyp$Df[2])
	}
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(terms(mod), "factors")
	intercept <- has.intercept(mod)
	p <- length(fixef(mod))
	I.p <- diag(p)
	assign <- attr(model.matrix(mod), "assign")
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <- names[-1]
	n.terms <- length(names)
	p <- teststat <- df <- rep(0, n.terms)
	for (i in 1:n.terms){
		hyp <- hyp.term(names[i])
		teststat[i] <- abs(hyp["statistic"])
		df[i] <- abs(hyp["df"])
		p[i] <- pchisq(teststat[i], df[i], lower.tail=FALSE) 
	}    
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c("Chisq", "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.III.lme <- function(mod, vcov., singular.ok=FALSE, ...){
	intercept <- has.intercept(mod)
	p <- length(coefficients(mod))
	I.p <- diag(p)
	names <- term.names(mod)
	n.terms <- length(names)
	assign <- attr(model.matrix(mod), "assign")
	df <- rep(0, n.terms)
	if (intercept) df[1] <- 1
	p <- teststat <-rep(0, n.terms)
	not.aliased <- !is.na(fixef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)        
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
		if (nrow(hyp.matrix) == 0){
			teststat[term] <- NA
			df[term] <- 0
			p[term] <- NA
		}
		else {
			hyp <- linearHypothesis(mod, hyp.matrix, 
					vcov.=vcov., singular.ok=singular.ok, ...)
			teststat[term] <-  hyp$Chisq[2] 
			df[term] <- abs(hyp$Df[2])
			p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
		}
	}
	result <- data.frame(teststat, df, p)
	row.names(result) <- names
	names(result) <- c ("Chisq",  "Df", "Pr(>Chisq)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

Anova.svyglm <- function(mod, ...) Anova.default(mod, ...)

Anova.rlm <- function(mod, ...) Anova.default(mod, test.statistic="F", ...)

Anova.coxme <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Wald", "LR"), ...){
    type <- as.character(type)
    type <- match.arg(type)
    test.statistic <- match.arg(test.statistic)
    switch(type,
        II=switch(test.statistic,
            LR=Anova.II.LR.coxme(mod, ...),
            Wald=Anova.default(mod, type="II", test.statistic="Chisq", ...)),
        III=switch(test.statistic,
            LR=stop("type-III LR tests not available for coxme models"),
            Wald=Anova.default(mod, type="III", test.statistic="Chisq", ...)),
        "2"=switch(test.statistic,
            LR=Anova.II.LR.coxme(mod, ...),
            Wald=Anova.default(mod, type="II", test.statistic="Chisq", ...)),
        "3"=switch(test.statistic,
            LR=stop("type-III LR tests not available for coxme models"),
            Wald=Anova.default(mod, type="III", test.statistic="Chisq")))
}

Anova.II.LR.coxme <- function(mod, ...){
    if (!requireNamespace("coxme")) stop("coxme package is missing")
    which.nms <- function(name) which(asgn == which(names == name))
    fac <-attr(terms(mod), "factors")
    names <- term.names(mod)
    n.terms <- length(names)
    if (n.terms < 2) return(anova(mod, test="Chisq"))
    X <- model.matrix(mod)
    asgn <- attr(X, 'assign')
    p <- LR <- rep(0, n.terms)
    df <- df.terms(mod)
    random <- mod$formulaList$random
    random <- sapply(random, as.character)[2, ]
    random <- paste(paste0("(", random, ")"), collapse=" + ")
    fixed <- as.character(mod$formulaList$fixed)[3]
    for (term in 1:n.terms){
        rels <- names[relatives(names[term], names, fac)]
        formula <- paste0(". ~ . - ", paste(c(names[term], rels), collapse=" - "), " + ", random)
        mod.1 <- update(mod, as.formula(formula))
        loglik.1 <- logLik(mod.1, type="integrated")
        mod.2 <- if (length(rels) == 0) mod
        else {
            formula <- paste0(". ~ . - ", paste(rels, collapse=" - "), " + ", random)
            update(mod, as.formula(formula))
        }
        loglik.2 <- logLik(mod.2, type="integrated")
        LR[term] <- -2*(loglik.1 - loglik.2)
        p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
    }
    result <- data.frame(LR, df, p)
    row.names(result) <- names
    names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
    result
}

