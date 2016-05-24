logbin <- function (formula, mono = NULL, data, subset, na.action, start = NULL, offset,
						control = list(...), model = TRUE, 
                        warn = TRUE, ...)
{
	call <- match.call()
	method <- "nplbin"
	
	family <- binomial(link = log)
	mu.eta <- family$mu.eta
	linkinv <- family$linkinv
	
	if(missing(data)) data <- environment(formula)
	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	control <- do.call("logbin.control", control)
	control2 <- control
	control2$trace <- (control$trace > 1)
	mt <- attr(mf, "terms")
	
	if (is.empty.model(mt)) stop("empty model")
	if (attr(mt, "intercept") != 1) stop("models without intercept are not supported by logbin")
	if (any(attr(mt, "order") > 1)) stop("models with interactions are not supported by logbin")
	if (attr(mt, "response") == 0) stop("missing response")
	
	Y <- model.response(mf, "numeric")
	if(length(dim(Y)) == 1L) {
		nm <- rownames(Y)
		dim(Y) <- NULL
		if(!is.null(nm)) names(Y) <- nm
	}
	
	allref <- logbin.allref(mt, mf, mono, start)
	design.numref <- sapply(allref$allref, length)
	
	offset <- as.vector(model.offset(mf))
	if(!is.null(offset)) {
		if(length(offset) != NROW(Y))
			stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
	
	best.model <- NULL
	best.loglik <- -Inf
	best.param <- NULL
	allconv <- TRUE
	
	if(length(allref$allref) == 0) {
		if(control$trace > 0) cat("logbin parameterisation 1/1\n")
		X <- model.matrix(allref$terms, allref$data)
		best.model <- nplbin(Y, X, offset, allref$start.new, control2)
		best.loglik <- best.model$loglik
		best.param <- 0
		allconv <- best.model$converged
		if(control$trace > 0 & control$trace <= 1)
			cat("Deviance =", best.model$deviance, "Iterations -", best.model$iter, "\n")
	} else {
		design.all <- expand.grid(lapply(design.numref, seq_len))
		nparam <- nrow(design.all)
		
		for(param in seq_len(nparam)) {
			if(control$trace > 0) cat("logbin parameterisation ",param,"/",nparam,"\n",sep="")
			X <- logbin.design(allref$terms, allref$data, allref$allref, design.all[param,])
			thismodel <- nplbin(Y, X, offset, if (param == 1) allref$start.new else NULL, control2)
			if(!thismodel$converged) allconv <- FALSE
			if(control$trace > 0 & control$trace <= 1)
				cat("Deviance =", thismodel$deviance, "Iterations -", thismodel$iter, "\n")
			if(thismodel$loglik > best.loglik) {
				best.model <- thismodel
				best.loglik <- thismodel$loglik
				best.param <- param
				if(thismodel$converged & !thismodel$boundary) break
			}
		}
	}
	
	if(warn) {
		if(!best.model$converged | (!allconv & best.model$boundary))
			warning(gettextf("%s: algorithm did not converge within %d iterations -- increase 'maxit'.",
                                method,control$maxit), 
                        call. = FALSE)
		if(best.model$boundary) {
			if(best.model$coefficients[1] > -control$bound.tol)
				warning(gettextf("%s: fitted probabilities numerically 1 occurred",
									method), call. = FALSE)
			else warning(gettextf("%s: MLE on boundary of parameter space",
									method), call. = FALSE)
		}
	}
	
	if(length(allref$allref) == 0) {
		np.coefs <- coefs <- best.model$coefficients
		nn.design <- design <- model.matrix(allref$terms, allref$data)
	} else {
		np.coefs <- best.model$coefficients
		nn.design <- logbin.design(allref$terms, allref$data, allref$allref, design.all[best.param,])
		reparam <- logbin.reparameterise(np.coefs, mt, mf, allref$allref, design.all[best.param,])
		coefs <- reparam$coefs
		design <- reparam$design
	}
	
	fit <- list(coefficients = coefs, residuals = best.model$residuals,
				fitted.values = best.model$fitted.values, rank = best.model$rank,
				family = family, linear.predictors = best.model$linear.predictors,
				deviance = best.model$deviance, loglik = best.model$loglik,
				aic = best.model$aic, aic.c = best.model$aic.c,
				null.deviance = best.model$null.deviance, iter = best.model$iter,
				prior.weights = best.model$prior.weights, weights = rep(1,NROW(Y)),
				df.residual = best.model$df.residual, df.null = best.model$df.null,
                y = best.model$y, x = design)
    if(model) fit$model <- mf
    fit2 <- list(converged = best.model$converged, boundary = best.model$boundary, 
                    na.action = attr(mf, "na.action"), call = call, formula = formula, 
                    terms = mt, data = data, offset = offset, control = control, 
                    method = method, xlevels = .getXlevels(mt, mf), 
                    xminmax = .getXminmax(mt, mf), np.coefficients = np.coefs,
                    nn.x = nn.design)
    
    fit <- c(fit, fit2)
	class(fit) <- c("logbin","glm","lm")
	
	fit
}