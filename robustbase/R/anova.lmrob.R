anova.lmrob <- function(object, ..., test = c("Wald", "Deviance"), verbose=getOption("verbose"))
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
	logical(length(dotargs))# FALSE
    else (names(dotargs) != "")
    if (any(named))
	warning("the following arguments to 'anova.lmrob' are invalid and \n",
		"dropped: ",
		pasteK(deparse(dotargs[named])))
    dotargs <- dotargs[!named]
    test <- match.arg(test)
    ## method argument has to end with 'M' (req. for refitting)
    if (test == "Deviance" && !grepl('M$', object$control$method))
      stop("For test = 'Deviance', the estimator chain has to end with 'M'")
    if (length(dotargs) > 0) {
	length.tl <- function(x) length(attr(terms(x),"term.labels"))
	isFormula <- vapply(dotargs, inherits, NA, what = "formula")
	h <- vapply(dotargs, length.tl, 0L)
	if(all(isFormula)) {
	    if(any(h >= length.tl(object)))
		stop("The first object does not contain the largest model")
	    modform <- dotargs
	}
	else {
	    if(verbose) message("All models are refitted except the largest one")
	    if(any(h >= length.tl(object))) {
		h <- c(length.tl(object),h)
		dotargs <- c(list(object), dotargs)[order(h, decreasing = TRUE)]
		object <- dotargs[[1]]
		if(!inherits(object, "lmrob"))
		    stop("anova.lmrob() only works for 'lmrob' objects")
		dotargs <- dotargs[-1]
	    }
	    modform <- lapply(dotargs, formula)
	}
	initCoef <- lapply(dotargs, coef)
	return(anovaLmrobList(object, modform, initCoef, test = test))
    }
    ##
    ## "'Anova Table' for a single model object
    stop("'Anova Table' for a single model not yet implemented")
}

anovaLmrobList <- function (object, modform, initCoef, test)
{
    responses <- as.character(lapply(modform, function(x) deparse(x[[2]])))
    if (!all(responses == deparse(formula(object)[[2]])))
	stop("Not the same response used in the fitted models")
    ##
    nobs <- length(object$residuals)
    nmodels <- length(modform) + 1
    tbl <- matrix(rep(NA, nmodels*4), ncol = 4)
    tbl[1,1] <- nobs[1] - length(coef(object))
    obj0 <- object
    for(k in 2:nmodels) {
	obj0 <- anovaLmrobPair(obj0, modform[[k-1]], initCoef[[k-1]],
			       test = test)
	tbl[k,] <- obj0$anova
	obj0$scale <- object$scale
    }
    ## return
    dimnames(tbl) <- list(1:nmodels,
			  c("pseudoDf", "Test.Stat", "Df", "Pr(>chisq)"))
    title <- switch(test,
		    Wald = "Robust Wald Test Table",
		    Deviance = "Robust Deviance Table",
		    stop("invalid 'test'"))
    variables <- c(list(formula(terms(object))), modform)
    topnote <- paste("Model ", format(1:nmodels), ": ", variables,
		     sep = "", collapse = "\n")
    note <- paste("Largest model fitted by lmrob(), i.e.",
                  object$control$method)
    ## paste("Models fitted by method '", methods[1], "'", sep="")
    structure(as.data.frame(tbl), heading = c(title, "", topnote, note,""),
	      class = c("anova", "data.frame"))
}

anovaLmrobPair <- function(FMfit, reduced.model, initCoef, test)
{
    ## 'FM': full model;  'RM' : reduced model
    X <- model.matrix(FMfit, data = FMfit$model)
    FMod <- FMfit$qr$pivot[1:FMfit$rank]
    asgn <- attr(X, "assign")
    FMt <- terms(FMfit)
    RMt <- terms(reduced.model)
    FMtl <- attr(FMt, "term.labels")
    RMtl <- attr(RMt, "term.labels")
    RMnumtl <- match(RMtl , FMtl, nomatch = -1)
    if(attr(RMt, "intercept") == 1) RMnumtl <- c(0, RMnumtl)
    if(any(is.na(match(RMnumtl, unique(asgn)))))
	stop("Models are not nested!")
    RMod0 <- seq(along = asgn)[!is.na(match(asgn, RMnumtl))]
    RMod <- intersect(RMod0, FMod)
    if (length(FMod) == length(RMod))
	stop("Models are not strictly nested")
    H0ind <- which(!FMod %in% RMod)
    H0coef <- coef(FMfit)[H0ind]
    df <- length(H0coef)
    pp <- FMfit$rank
    switch(test, "Wald" = {
	t.cov <- FMfit$cov
	t.chisq <- sum(H0coef * solve(t.cov[H0ind, H0ind], H0coef))
        ## return
	c(FMfit,
          list(anova = c(nrow(X)-pp+df, t.chisq, df,
               pchisq(as.vector(t.chisq), df = df, lower.tail = FALSE))))
    },
    "Deviance" = {
	y <- FMfit$residuals + FMfit$fitted.values
	s0 <- FMfit$scale
	psi <- function(u, deriv = 0)
	    Mpsi(u, cc = FMfit$control$tuning.psi,
                   psi = FMfit$control$psi, deriv)
	iC <-
	    if(is.null(initCoef)) {
		res <- as.vector(y - X[,RMod] %*% FMfit$coef[RMod])
		psiRes <- psi(res/s0)
		if(sum(abs(psiRes) < 1e-08) > 0.6*nrow(X))
		    stop("Please fit the nested models by lmrob")
		FMfit$coef[RMod]
	    } else {
                idx <- !is.na(initCoef)
                if (any(idx != RMod0 %in% RMod))
                    stop("NA coefs in full and reduced model do not match")
                initCoef[idx]
            }

	RMfit <- lmrob..M..fit(x = X[,RMod, drop=FALSE], y = y,
			       beta.initial = iC, scale = s0,
			       control = FMfit$control)
	FMres <- as.vector(y - X[,FMod] %*% FMfit$coef[FMod])
	RMres <- RMfit$resid ## as.vector(y - X[,RMod] %*% RMfit$coef)
	FM_sRho <- sum(psi(FMres/s0, deriv = -1))
	RM_sRho <- sum(psi(RMres/s0, deriv = -1))
	tauStar <- mean(psi(FMres/s0,	deriv = 1)) /
		   mean(psi(FMres/s0)^2, deriv = 0)
	t.chisq <- 2*tauStar*(RM_sRho - FM_sRho)
	## return
	c(RMfit,
          list(anova = c(nrow(X)-pp+df, t.chisq, df,
               pchisq(as.vector(t.chisq), df = df, lower.tail = FALSE))))
    },
    stop("test ", test, " not yet implemented"))
} ## anovaLmrobPair

