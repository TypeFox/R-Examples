addreg <- function (formula, mono = NULL, family, data, standard, subset, na.action, start = NULL, offset,
                        control = list(...), model = TRUE, warn = TRUE, ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
		if (identical(family, negbin1)) family <- family(link = "identity", phi = NA)
        else family <- family(link = "identity")
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
    if (family$link!="identity" | !(family$family %in% c("poisson","binomial") | substr(family$family,1,7) == "negbin1"))
        stop("family(link) must be one of: poisson(identity), binomial(identity), negbin1(identity)")
    
    if (family$family == "poisson") method <- "nnpois"
    else if (substr(family$family,1,7) == "negbin1") method <- "nnnegbin"
    else if (family$family == "binomial") method <- "addbin"
    
    if (missing(data)) data <- environment(formula)
    
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "offset", "standard"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    control <- do.call("addreg.control", control)
    control2 <- control
    control2$trace <- (control$trace > 1)
    mt <- attr(mf, "terms")
    
    if (is.empty.model(mt)) stop("empty model")
    if (attr(mt,"intercept") != 1) stop("models without intercept are not supported by addreg")
    if (any(attr(mt,"order") > 1)) stop("models with interactions are not supported by addreg")
    if (attr(mt,"response") == 0) stop("missing response")
    
    Y <- model.response(mf, "numeric")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
    }
       
    allref <- addreg.allref(mt, mf, mono, family, start)
	design.numref <- sapply(allref$allref, length)
    
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of values in 'offset' is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
        if (family$family == "binomial")
            warning("'offset' is not supported for binomial family", call. = FALSE)
    }
    standard <- as.vector(model.extract(mf,"standard"))
    if(!is.null(standard)) {
        if (length(standard) != NROW(Y))
            stop(gettextf("number of values in 'standard' is %d should equal %d (number of observations)",
                length(standard), NROW(Y)), domain = NA)
        if (family$family == "binomial")
            warning("'standard' is not used for binomial family", call. = FALSE)
		if (any(standard <= 0))
			stop("standard must be positive")
	}
    best.model <- NULL
    best.loglik <- -Inf
    best.param <- NULL
    allconv <- TRUE
    
	if (length(allref$allref) == 0) {
        if (control$trace > 0) cat(method,"parameterisation 1/1\n")
		X <- model.matrix(allref$terms, allref$data)
        if (family$family == "poisson")
            best.model <- nnpois(Y, X, standard, offset, allref$start.new, control2)
        else if (substr(family$family,1,7) == "negbin1")
            best.model <- nnnegbin(Y, X, standard, offset, allref$start.new, control2)
        else if (family$family == "binomial")
            best.model <- addbin(Y, X, allref$start.new, control, allref)
        best.loglik <- best.model$loglik
        best.param <- 0
        allconv <- best.model$converged
        if (control$trace > 0 & control$trace <= 1)
			if (substr(family$family,1,7) == "negbin1") cat("Log-likelihood =",best.model$loglik,
                                                            "Iterations -",best.model$iter,"\n")
			else if (method != "addbin")
                cat("Deviance =", best.model$deviance, "Iterations -", best.model$iter, "\n")
	} else {
		design.all <- expand.grid(lapply(design.numref, seq_len))
		nparam <- nrow(design.all)
        
		for(param in seq_len(nparam)) {
            if (control$trace > 0) cat(method," parameterisation ",param,"/",nparam,"\n",sep="")
			X <- addreg.design(allref$terms, allref$data, allref$allref, design.all[param,])
            if (family$family == "poisson")
                thismodel <- nnpois(Y, X, standard, offset, if (param == 1) allref$start.new else NULL, control2)
            else if (substr(family$family,1,7) == "negbin1")
                thismodel <- nnnegbin(Y, X, standard, offset, if (param == 1) allref$start.new else NULL, control2)
            else if (family$family == "binomial")
                thismodel <- addbin(Y, X, if (param == 1) allref$start.new else NULL, control, allref)
            if (!thismodel$converged) allconv <- FALSE
            if (control$trace > 0 & control$trace <= 1) 
				if (substr(family$family,1,7) == "negbin1") cat("Log-likelihood =",thismodel$loglik,
                                                                "Iterations -",thismodel$iter,"\n")
				else if (method != "addbin")
                    cat("Deviance =",thismodel$deviance,"Iterations -",thismodel$iter,"\n")
            if (thismodel$loglik > best.loglik) {
                best.model <- thismodel
                best.loglik <- thismodel$loglik
                best.param <- param
                if(thismodel$converged & !thismodel$boundary) break
            }
		}
	}
    
    if (warn) {
        if (!best.model$converged | (!allconv & best.model$boundary)) 
            warning(gettextf("%s: algorithm did not converge within %d iterations -- increase 'maxit'.",
                                method,control$maxit), 
                        call. = FALSE)
        if (best.model$boundary) {
            if (best.model$coefficients[1] < control$bound.tol) {
                if (family$family == "poisson" || substr(family$family,1,7) == "negbin1")
                    warning(gettextf("%s: fitted rates numerically 0 occurred",
                                        method), call. = FALSE)
                else if (family$family == "binomial")
                    warning(gettextf("%s: fitted probabilities numerically 0 or 1 occurred",
                                        method), call. = FALSE)
            } else warning(gettextf("%s: MLE on boundary of parameter space",
                                        method), call. = FALSE)
        }  
    }
    
    if (length(allref$allref) == 0) {
        nn.coefs <- coefs <- best.model$coefficients
        nn.design <- design <- model.matrix(allref$terms, allref$data)
    } else {
        nn.coefs <- best.model$coefficients
        nn.design <- addreg.design(allref$terms, allref$data, allref$allref, design.all[best.param,])
		reparam <- addreg.reparameterise(nn.coefs, mt, mf, allref$allref, design.all[best.param,])
        coefs <- reparam$coefs
        design <- reparam$design
    }
    
    fit <- list(coefficients = coefs)
    if (substr(family$family,1,7) == "negbin1") fit$scale <- best.model$scale
    
    fit2 <- list(residuals = best.model$residuals,
                fitted.values = best.model$fitted.values, 
                rank = best.model$rank, family = best.model$family, 
                linear.predictors = best.model$linear.predictors,
                deviance = best.model$deviance, loglik = best.model$loglik,
                aic = best.model$aic, aic.c = best.model$aic.c,
                null.deviance = best.model$null.deviance, iter = best.model$iter,
                prior.weights = best.model$prior.weights, weights = rep(1,NROW(Y)),
                df.residual = best.model$df.residual, df.null = best.model$df.null,
                y = best.model$y, x = design)
    if (model) {
        fit2$model <- mf
        if (family$family == "binomial") fit2$model.addpois <- best.model$model.addpois
    }
    fit3 <- list(converged = best.model$converged, boundary = best.model$boundary, 
                na.action = attr(mf, "na.action"), call = call, formula = formula, 
                terms = mt, data = data, offset = best.model$offset, standard = best.model$standard, 
                control = control, method = method, xlevels = .getXlevels(mt, mf),
                xminmax = .getXminmax(mt, mf))
    if (family$family == "poisson" || substr(family$family,1,7) == "negbin1") {
        fit3$nn.coefficients = nn.coefs
        fit3$nn.x = nn.design
    }
    fit <- c(fit, fit2, fit3)
    class(fit) <- c("addreg","glm","lm")
    
    fit
}