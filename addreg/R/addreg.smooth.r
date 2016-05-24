addreg.smooth <- function (formula, mono = NULL, family, data, standard, subset, na.action,
					offset, control = list(...), model = TRUE, model.addreg = FALSE, ...) {
	call <- match.call()
	
	if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family))
        if (identical(family, negbin1)) {
            family <- family(link = "identity", phi = NA)
            famname <- "negbin1"
        }
        else {
            family <- family(link = "identity")
            famname <- family$family
        }
    if(is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
    if(family$link!="identity" | !(family$family %in% c("poisson","binomial") | substr(family$family,1,7) == "negbin1"))
        stop("family(link) must be one of: poisson(identity), binomial(identity), negbin1(identity)")
	
	if(missing(data)) data <- environment(formula)
	
	gp <- interpret.addreg.smooth(formula)
	mf <- match.call(expand.dots = FALSE)
	mf$formula <- gp$fake.formula
    m <- match(c("formula", "data", "subset", "offset", "standard"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mdata <- mf
    mdata[[1L]] <- as.name("get_all_vars")
    mdata <- eval(mdata, parent.frame())
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
	if (!missing(na.action)) mf$na.action <- na.action
    mf <- eval(mf, parent.frame())
	control <- do.call("addreg.control", control)
	control2 <- control
	control2$trace <- pmax(0, as.numeric(control$trace) - 1)
	mt <- attr(mf, "terms")
    
    os <- model.offset(mf)
    std <- model.extract(mf, "standard")
    
    allknots <- expand.grid(lapply(gp$smooth.spec,"[[","knot.range"))
    n.allknots <- nrow(allknots)
    
	bestk <- NULL
    bestk.model <- NULL
    bestk.aic <- Inf
	bestk.allref <- NULL
    bestk.param <- NULL
	bestk.knots <- NULL
    
    for(k in seq_len(n.allknots)) {
        if(control$trace > 0)
            cat("Knots: ",paste(allknots[k,],collapse=", "),"\n",sep="")
        allref <- addreg.smooth.allref(mt, mdata, mono, family, gp, allknots[k,])
        design.numref <- sapply(allref$allref, length)
        design.all <- expand.grid(lapply(design.numref, seq_len))
        nparam <- nrow(design.all)
        
        best.model <- NULL
        best.loglik <- -Inf
        best.param <- NULL
		best.knots <- NULL
        allconv <- TRUE
        for (param in seq_len(nparam)) {
            if(control$trace > 1) cat("addreg.smooth parameterisation ",param,"/",nparam,"\n",sep="")
            modelspec <- addreg.smooth.design(gp, allref, allknots[k,,drop=FALSE], 
                            design.all[param,,drop=FALSE])
            data.new <- modelspec$data
			data.new[["(offset)"]] = os
			data.new[["(standard)"]] = std
			modelf <- call("addreg",formula = eval(modelspec$formula), mono = eval(modelspec$monotonic), 
                        family = famname, data = as.name("data.new"), control = control2, 
						warn = FALSE, fit = TRUE)
			if (!is.null(os)) modelf$offset <- as.name("(offset)")
			if (!is.null(std)) modelf$standard <- as.name("(standard)")
			if (!missing(subset)) modelf$subset <- subset
			if (!missing(na.action)) modelf$na.action <- na.action
            modelf$model <- TRUE
			thismodel <- eval(modelf)
            if(!thismodel$converged) allconv <- FALSE
			if(thismodel$loglik > best.loglik) {
				best.model <- thismodel
				best.loglik <- thismodel$loglik
				best.param <- design.all[param,,drop=FALSE]
				best.knots <- modelspec$knots
				if(thismodel$converged & !thismodel$boundary) break
			}
        }
		if(!best.model$converged || (!allconv & best.model$boundary))
			warning(gettextf("%s: algorithm did not converge within %d iterations -- increase 'maxit'.",
							best.model$method,control$maxit), 
					call. = FALSE)
		if(control$trace > 0)
			cat("AIC_c:",best.model$aic.c,"\n")
		if(best.model$aic.c < bestk.aic) {
			bestk <- k
			bestk.model <- best.model
			bestk.aic <- best.model$aic.c
			bestk.allref <- allref
			bestk.param <- best.param
			bestk.knots <- best.knots
		}
    }
	
	if(bestk.model$boundary) {
		if(family$family == "poisson" || substr(family$family,1,7) == "negbin1") {
			if(bestk.model$nn.coefficients[1] < control$bound.tol)
				warning(gettextf("%s: fitted rates numerically 0 occurred",
									bestk.model$method), call. = FALSE)
			else
				warning(gettextf("%s: MLE on boundary of parameter space",
									bestk.model$method), call. = FALSE)
		} else if(family$family == "binomial") {
			if(bestk.model$model.addpois$nn.coefficients[1] < control$bound.tol)
				warning(gettextf("%s: fitted probabilities numerically 0 or 1 occurred",
									bestk.model$method), call. = FALSE)
			else
				warning(gettextf("%s: MLE on boundary of parameter space",
									bestk.model$method), call. = FALSE)
		}
	}  
	
	reparam.call <- call("addreg.smooth.reparameterise", coefficients = bestk.model$coefficients,
							interpret = gp, allref = bestk.allref, knots = bestk.knots, 
							design.knots = allknots[bestk,,drop=FALSE], design.param = bestk.param)

	if(!missing(subset)) reparam.call$subset <- subset
	if(!missing(na.action)) reparam.call$na.action <- na.action
	reparam <- eval(reparam.call)
	
	fit <- list(coefficients = reparam$coefs)
	if(substr(family$family,1,7) == "negbin1") fit$scale <- bestk.model$scale
	
	fit2 <- list(residuals = bestk.model$residuals,
				fitted.values = bestk.model$fitted.values,
				rank = bestk.model$rank, family = bestk.model$family,
				linear.predictors = bestk.model$linear.predictors,
				deviance = bestk.model$deviance, loglik = bestk.model$loglik,
				aic = bestk.model$aic, aic.c = bestk.model$aic.c,
				null.deviance = bestk.model$null.deviance, iter = bestk.model$iter,
				prior.weights = bestk.model$prior.weights, weights = bestk.model$weights,
				df.residual = bestk.model$df.residual, df.null = bestk.model$df.null,
                y = bestk.model$y, x = reparam$design)
    if(model) fit2$model <- reparam$mf
    if(model.addreg) fit2$model.addreg <- bestk.model
    xminmax.smooth <- bestk.model$xminmax
    xminmax.smooth[reparam$smoothnames] <- NULL
	fit3 <- list(converged = bestk.model$converged, boundary = bestk.model$boundary,
				na.action = attr(reparam$mf, "na.action"), call = call, formula = formula,
				full.formula = gp$full.formula,
                terms = mt, terms.full = reparam$mt, data = data, offset = os, standard = bestk.model$standard, 
                control = control, method = bestk.model$method, xlevels = bestk.model$xlevels,
                xminmax = xminmax.smooth, knots = bestk.knots)
                
	fit <- c(fit, fit2, fit3)
	class(fit) <- c("addreg.smooth","addreg","glm","lm")
    
    fit
}