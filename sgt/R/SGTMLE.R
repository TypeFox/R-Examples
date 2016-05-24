print.sgtest = function(x, ...) {
	cat("Skewed Generalized T MLE Fit\n")
	cat("Best Result with", x$best.method.used, "Maximization\n")
	cat("Convergence Code ", x$convcode, ": ", sep="")
	if (x$convcode == 0) {
		cat("Successful Convergence\n")
	} else if (x$convcode == 1) {
		cat("Iteration Limit Reached\n")
	} else if (x$convcode == 20) {
		cat("Initial Set of Parameters is not Permissible\n")
	} else {
		cat("Maximization Failure\n")
	}
	cat("Iterations:", x$niter)
	cat(", Log-Likelihood:", x$maximum, "\n")
	cat("Estimate(s):\n")
	print(noquote(format(round(x$estimate, 4), nsmall=4)))
}

summary.sgtest = function(object, ...) {
	object$z.score = object$estimate/object$std.error
	object$p.value = 2*stats::pnorm(-abs(object$z.score))
	names(object$z.score) = names(object$estimate)
	names(object$p.value) = names(object$estimate)
	object$summary.table = as.data.frame(cbind(object$estimate, object$std.error, object$z.score, object$p.value), row.names=names(object$estimate))
	names(object$summary.table) = c("Est.", "Std. Err.", "z", "P>|z|")
	class(object) = c("summary.sgtest", "list")
	object
}

print.summary.sgtest = function(x, ...) {
	cat("Skewed Generalized T MLE Fit\n")
	cat("Best Result with", x$best.method.used, "Maximization\n")
	cat("Convergence Code ", x$convcode, ": ", sep="")
	if (x$convcode == 0) {
		cat("Successful Convergence\n")
	} else if (x$convcode == 1) {
		cat("Iteration Limit Reached\n")
	} else if (x$convcode == 20) {
		cat("Initial Set of Parameters is not Permissible\n")
	} else {
		cat("Maximization Failure\n")
	}
	cat("Iterations:", x$niter)
	cat(", Log-Likelihood:", x$maximum, "\n\n")
	signif = rep("", length.out = length(x$estimate))
	signif[x$p.value < 0.1] = "."
	signif[x$p.value < 0.05] = "*"
	signif[x$p.value < 0.01] = "**"
	signif[x$p.value < 0.001] = "***"
	outtable = cbind(format(round(x$summary.table, 4), nsmall = 4), signif)
	names(outtable)[5] = ""
	print(noquote(outtable))
	cat("---\n")
	cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
}

sgt.mle = function(X.f, mu.f = mu ~ mu, sigma.f = sigma ~ sigma, lambda.f = lambda ~ lambda, p.f = p ~ p, q.f = q ~ q, data = parent.frame(), start, subset, method = c("Nelder-Mead", "BFGS"), itnmax = NULL, hessian.method="Richardson", gradient.method="Richardson", mean.cent = TRUE, var.adj = TRUE, ...) {
	formList = list(X = X.f, mu = mu.f, sigma = sigma.f, lambda = lambda.f, p = p.f, q = q.f)
	varNames = NULL
	envir = new.env()
	for(i in 1:6) {
		formList[[i]] = stats::as.formula(formList[[i]])
		if(length(formList[[i]]) == 2L) {
			formList[[i]][[3L]] = formList[[i]][[2L]]
			formList[[i]][[2L]] = as.name(names(formList)[i])
    	} else if(as.character(formList[[i]][[2L]]) != names(formList)[i]) {
    		warning(paste('The left hand side of ', names(formList)[i], '.f was changed from ',as.character(formList[[i]][[2L]]),' to ', names(formList)[i],sep=''))
    	}
    	varNames = c(varNames,all.vars(formList[[i]][[3L]]))
	}
	if (class(data)[1L] == "matrix") data = as.data.frame(data)
	if (!is.list(data) && !is.environment(data)) 
        stop("'data' must be a list or an environment")
    start = as.list(start)
    if(is.null(names(start)))
    	stop("'start' must be a named list or named numeric vector")
    if("" %in% names(start))
    	stop("at least one of the elements in 'start' is missing a name")
	parNames = names(start)
	varNames = varNames[is.na(match(varNames,parNames))]
	if(length(varNames) == 0L)
		stop("there is no reference to data in the given formulas")
	for(i in varNames) {
		if(!exists(i, data))
			stop(paste(i,"is not contained in 'start' and it is not found in 'data'"))
		assign(i, eval(parse(text=paste("as.numeric(data$",i,")",sep=""))), envir)
	}
	if(length(varNames) > 1) {
		for(i in 2:length(varNames)) {
			if(length(eval(parse(text=paste("envir$", varNames[1L], sep="")))) != length(eval(parse(text=paste("envir$", varNames[i], sep="")))))
				stop(paste("the length of the variable", varNames[i], "does not match the length of the variable", varNames[1L]))
		}
	}
	control = list(...)
	if(!is.null(control$maximize)) stop("'maximize' option not allowed")
	if(!missing(subset))
		for(i in varNames) assign(i, eval(parse(text=paste("envir$", i, "[subset]", sep=""))), envir)
	keep = rep(TRUE,length(eval(parse(text=paste("envir$", varNames[1L], sep="")))))
	for(i in varNames) keep = keep & is.finite(eval(parse(text=paste("envir$", i, sep=""))))
	for(i in varNames) assign(i, eval(parse(text=paste("envir$", i, "[keep]", sep=""))), envir)
	loglik = function(params) {
		for (i in 1:length(parNames)) assign(parNames[i], unlist(params[i]))
        X = eval(formList[[1L]][[3L]])
        mu = eval(formList[[2L]][[3L]])
        sigma = eval(formList[[3L]][[3L]])
        lambda = eval(formList[[4L]][[3L]])
        p = eval(formList[[5L]][[3L]])
        q = eval(formList[[6L]][[3L]])
		sum(dsgt(X, mu, sigma, lambda, p, q, mean.cent, var.adj, log=TRUE))
	}
	environment(loglik) = envir
	negloglik = function(params) {-loglik(params)} 
	if (!is.finite(loglik(start))) stop("'start' yields infinite or non-computable SGT function values")
	optimum = suppressWarnings(optimx::optimx(par = unlist(start), fn = negloglik, method = method, itnmax = itnmax, control = control))
	minimum = min(optimum$value, na.rm = TRUE)
	if(!is.finite(minimum)) stop("All Maximization Methods Failed")
	whichbest = max(which(minimum == optimum$value))
	optimal = optimum[whichbest,]
	estimate = as.numeric(optimum[whichbest, 1:length(parNames)])
	names(estimate) = parNames
	H = tryCatch(numDeriv::hessian(loglik, estimate, method = hessian.method), error = function(e) {warning("hessian matrix calculation failed"); return(as.matrix(NaN))})
	varcov = tryCatch(-qr.solve(H), error = function(e) {warning("covariance matrix calculation failed due to a problem with the hessian"); return(as.matrix(NaN))})
	std.error = sqrt(diag(varcov))
	if(is.finite(varcov[1,1])) names(std.error) = parNames
	gradient = tryCatch(numDeriv::grad(loglik, estimate, method = gradient.method), error = function(e) {warning("gradient calculation failed"); return(NaN)})
	result = list(maximum = -minimum, estimate = estimate, convcode = as.numeric(optimal$convcode), niter = as.numeric(optimal$niter), best.method.used = row.names(optimal), optimx = optimum, hessian = H, gradient = gradient, varcov = varcov, std.error = std.error)
	class(result) = c("sgtest", class(result))
	return(result)
}

