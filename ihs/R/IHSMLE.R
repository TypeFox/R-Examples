summary.MLE = function(object, ...) {
	params = object$parameters
	class(object) = c("maxLik","maxim","list")
	object = summary(object)
	object$parameters = params
	class(object) = c("summary.MLE", "summary.maxLik")
	object
}

print.summary.MLE = function(x, ...) {
	cat('--------------------------------------------\n')
	cat('Maximum Likelihood estimation\n')
	cat(paste(x$maximType, ", ", x$iterations, " iterations\n", sep=""))
	cat(paste("Return code ", x$returnCode, ": ", x$returnMessage, "\n", sep=""))
	cat(paste("Log-Likelihood: ", format(x$loglik, digits=4, nsmall=4), "\n",sep=""))
	cat(paste(x$NActivePar,"free parameter"))
	if(x$NActivePar > 1) cat("s\n")
	cat("Estimates:\n")
	Parameters = x$parameters
	print(cbind(Parameters,format(as.data.frame(x$estimate), digits=4)))
	cat('--------------------------------------------\n')
}

print.mult.MLE = function(x, ...) {
	cat('Maximum Likelihood estimation\n')
	maxval = x[[1L]]$maximum
	maxindex = 1L
	for (i in 2:length(x)) {
		if (x[[i]]$maximum > maxval) {
			maxval = x[[i]]$maximum
			maxindex = i
		}
	}
	subx = x[[maxindex]]
	cat('Highest Log-Likelihood value: ')
	cat(paste(format(subx$maximum, digits=4, nsmall=4),"\n"))
	cat(paste(subx$type, ", ", subx$iterations, " iterations\n", sep=""))
	cat(paste("Return code ", subx$code, ": ", subx$message, "\n", sep=""))
	cat("Estimate(s): ")
	cat(format(subx$estimate, digits=4, nsmall=4))
	cat("\n")
}

summary.mult.MLE = function(object, ...) {
	maxNames = names(object)
	result = list()
	for (i in maxNames) {
		eval(parse(text=paste("result$", i, " = summary(object$", i, ")", sep="")))
	}
    class(result) = c("summary.mult.MLE", class(result))
    result
}

print.summary.mult.MLE = function(x, ...) {
    cat('Maximum Likelihood estimation\n')
    cat('--------------------------------------------\n')
    for (i in 1:length(x)) {
        cat(paste(x[[i]]$maximType, ", ", x[[i]]$iterations, " iterations\n", sep=""))
        cat(paste("Return code ", x[[i]]$returnCode, ": ", x[[i]]$returnMessage, "\n", sep=""))
        cat(paste("Log-Likelihood: ", format(x[[i]]$loglik, digits=4, nsmall=4), "\n",sep=""))
        cat(paste(x[[i]]$NActivePar,"free parameter"))
        if(x[[i]]$NActivePar > 1) cat("s\n")
            cat("Estimates:\n")
        Parameters = x[[i]]$parameters
        print(cbind(Parameters, format(as.data.frame(x[[i]]$estimate), digits=4)))
        cat('--------------------------------------------\n')
    }
}

ihs.mle = function(X.f, mu.f = mu ~ mu, sigma.f = sigma ~ sigma, lambda.f = lambda ~ lambda, k.f = k ~ k, data = parent.frame(), start, subset, method = 'BFGS', constraints = NULL, follow.on = FALSE, iterlim = 5000, ...) {
	formList = list(X = X.f, mu = mu.f, sigma = sigma.f, lambda = lambda.f, k = k.f)
	varNames = NULL
	envir = new.env()
	for(i in 1:5) {
		formList[[i]] = as.formula(formList[[i]])
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
				stop(paste("The length of the variable", varNames[i], "does not match the length of the variable", varNames[1L]))
		}
	}
	if(!missing(subset))
		for(i in varNames) assign(i, eval(parse(text=paste("envir$", i, "[subset]", sep=""))), envir)
	keep = rep(TRUE,length(eval(parse(text=paste("envir$", varNames[1L], sep="")))))
	for(i in varNames) keep = keep & is.finite(eval(parse(text=paste("envir$", i, sep=""))))
	for(i in varNames) assign(i, eval(parse(text=paste("envir$", i, "[keep]", sep=""))), envir)
	logLik = function(params) {
		for (i in 1:length(parNames)) assign(parNames[i], params[i])
        X = eval(formList[[1L]][[3L]])
        mu = eval(formList[[2L]][[3L]])
        sigma = eval(formList[[3L]][[3L]])
        lambda = eval(formList[[4L]][[3L]])
        k = eval(formList[[5L]][[3L]])
		dihs(X, mu, sigma, lambda, k, log=TRUE)
	}
	environment(logLik) = envir
	method = toupper(method)
	if (sum(is.na(match(method, c("NR", "BFGS", "BHHH", "SANN", "CG", "NM")))) > 0)
		stop("'method' argument is not valid")
	if (length(method) == 1) {
		result = maxLik(logLik = logLik, start = as.numeric(unlist(start)), method = method, constraints = constraints, iterlim = iterlim, ...)
		class(result) = c("MLE", class(result))
		result$parameters = parNames
	} else {
		result = list()
		envir2 = new.env()
		for (i in 1:length(method)) {
			if(!exists(method[i], envir2))
				assign(method[i], 1, envir2)
			else
				assign(method[i], eval(parse(text=paste("envir2$", method[i], sep="")))+1, envir2)
			resultStr = paste("result$", method[i], sep="")
			if (eval(parse(text=paste("envir2$", method[i], sep=""))) > 1)
				resultStr = paste(resultStr, eval(parse(text=paste("envir2$", method[i], sep=""))), sep="")
			if (length(iterlim) == 1)
				iters = iterlim
			else if (length(iterlim) == length(method))
				iters = iterlim[i]
			else stop("'iterlim' must either be a positive integer or a vector of positive integers with the same length as 'method'")
			out = tryCatch(maxLik(logLik = logLik, start = as.numeric(unlist(start)), method = method[i], constraints = constraints, iterlim = iters, ...), error = function(e) e)
			if (!("error" %in% class(out))) {
				class(out) = c("MLE", class(out))
				out$parameters = parNames
				eval(parse(text=paste(resultStr, " = out", sep="")))
				if (follow.on)
					start = out$estimate
			} else {
				warning(paste("maximisation using the ", method[i], " method failed:\n", out$message, "; this method must be skipped", sep=""))
			}
		}
		if (length(result) == 0)
			stop("all maximisation methods failed")
		class(result) = c("mult.MLE", class(result))
	}
	result
}