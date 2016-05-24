#############################################################
#                                                           #
#	anova.wle.glm.root function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: June, 16, 2009                                #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2009 Claudio Agostinelli              #
#                                                           #
#############################################################

## Function developed from 'anova.glm' R version 2.10.0
anova.wle.glm.root <- function(object, ..., dispersion=NULL, test=NULL) {
    warnings('Hey Claudio: remember to finish this function before to release the package!')
    if (class(object)!='wle.glm.root')
      stop("Use 'extractRoot(object)' to extract a single root from a wle.glm object")
    ## check for multiple objects
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
	rep(FALSE, length(dotargs)) else (names(dotargs) != "")
    if(any(named))
	warning("the following arguments to 'anova.wle.glm.root' are invalid and dropped: ", paste(deparse(dotargs[named]), collapse=", "))
    dotargs <- dotargs[!named]
    is.glm <- unlist(lapply(dotargs,function(x) inherits(x,"wle.glm.root")))
    dotargs <- dotargs[is.glm]
    nargs <- length(dotargs)
    if (nargs) 
      return(anova.wleglmlist(c(list(object), dotargs), dispersion = dispersion, test = test))
    
    ## extract variables from model
    varlist <- attr(object$terms, "variables")
    ## must avoid partial matching here.
    x <- if (n <- match("x", names(object), 0L))
	   object[[n]]
	 else model.matrix(object)
    varseq <- attr(x, "assign")
    nvars <- max(0, varseq)
    resdev <- resdf <- NULL

    ## if there is more than one explanatory variable then
    ## recall glm.fit to fit variables sequentially

    if(nvars > 1) {
	method <- object$method
	if(!is.function(method))
	    method <- get(method, mode = "function", envir=parent.frame())
        ## allow for 'y = FALSE' in the call (PR#13098)
        y <- object$y
        if(is.null(y)) { ## code from residuals.glm
            mu.eta <- object$family$mu.eta
            eta <- object$linear.predictors
            y <-   object$fitted.values + object$residuals * mu.eta(eta)
        }
        nobs <- NROW(y)
	for(i in 1L:(nvars-1)) {
	    ## explanatory variables up to i are kept in the model
	    ## use method from glm to find residual deviance
	    ## and df for each sequential fit
	    fit <- method(x=x[, varseq <= i, drop = FALSE],
			  y=y,
			  weights=object$prior.weights*object$wle.weights,
			  start	 =object$start,
			  offset =object$offset,
			  family =object$family,
			  control=object$control)
	    resdev <- c(resdev, fit$deviance)
	    resdf <- c(resdf, fit$df.residual + sum(object$wle.weights) - nobs)
	}
    }

    ## add values from null and full model

    resdf <- c(object$df.null, resdf, object$df.residual)
    resdev <- c(object$null.deviance, resdev, object$deviance)

    ## construct table and title

    table <- data.frame(c(NA, -diff(resdf)),
			c(NA, pmax(0, -diff(resdev))), resdf, resdev)
    tl <- attr(object$terms, "term.labels")
    if (length(tl) == 0L) table <- table[1,,drop=FALSE] # kludge for null model
    dimnames(table) <- list(c("NULL", tl),
			    c("Df", "Deviance", "Resid. Df", "Resid. Dev"))
    title <- paste("Robust Analysis of Deviance Table", "\n\nModel: ",
		   object$family$family, ", link: ", object$family$link,
		   "\n\nResponse: ", as.character(varlist[-1L])[1L],
		   "\n\nTerms added sequentially (first to last)\n\n", sep="")

    ## calculate test statistics if needed

    df.dispersion <- Inf
    if(is.null(dispersion)) {
	dispersion <- summary(object, dispersion=dispersion)$dispersion
	df.dispersion <- if (dispersion == 1) Inf else object$df.residual
    }
    if(!is.null(test)) {
        if(test == "F" && df.dispersion == Inf) {
            fam <- object$family$family
            if(fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a %s family is inappropriate",
                                 fam),
                        domain = NA)
            else
                warning("using F test with a fixed dispersion is inappropriate")
        }
	table <- stat.anova(table=table, test=test, scale=dispersion,
			    df.scale=df.dispersion, n=sum(object$wle.weights))
    }
    structure(table, heading = title, class= c("anova", "data.frame"))
}


anova.wleglmlist <- function(object, ..., dispersion=NULL, test=NULL) {

    ## find responses for all models and remove
    ## any models with a different response

    responses <- as.character(lapply(object, function(x) {
	deparse(formula(x)[[2L]])} ))
    sameresp <- responses==responses[1L]
    if(!all(sameresp)) {
	object <- object[sameresp]
	warning("models with response ", deparse(responses[!sameresp]),
		" removed because response differs from model 1")
    }

    ns <- sapply(object, function(x) length(x$residuals))
    if(any(ns != ns[1L]))
	stop("models were not all fitted to the same size of dataset")

    ## calculate the number of models

    nmodels <- length(object)
    if(nmodels==1)
	return(anova.wle.glm.root(object[[1L]], dispersion=dispersion, test=test))

    ## extract statistics

    resdf  <- as.numeric(lapply(object, function(x) x$df.residual))
    resdev <- as.numeric(lapply(object, function(x) x$deviance))

    ## construct table and title

    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)),
			c(NA, -diff(resdev)) )
    variables <- lapply(object, function(x)
			paste(deparse(formula(x)), collapse="\n") )
    dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", "Df",
					 "Deviance"))
    title <- "Robust Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1L:nmodels),": ",
		     variables, sep="", collapse="\n")

    ## calculate test statistic if needed

    if(!is.null(test)) {
        warning('Hey Claudio: the order to choose the bigmodel is wrong fix me!')
	bigmodel <- object[[order(resdf)[1L]]]
	dispersion <- summary(bigmodel, dispersion=dispersion)$dispersion
	df.dispersion <- if (dispersion == 1) Inf else min(resdf)
        if(test == "F" && df.dispersion == Inf) {
            fam <- bigmodel$family$family
            if(fam == "binomial" || fam == "poisson")
              warning(gettextf("using F test with a '%s' family is inappropriate", fam), domain = NA, call. = FALSE)
            else
              warning("using F test with a fixed dispersion is inappropriate")
        }
	table <- stat.anova(table = table, test = test, scale = dispersion, df.scale = df.dispersion, n = sum(bigmodel$wle.weights))
    }
    structure(table, heading = c(title, topnote), class = c("anova", "data.frame"))
}
