########################################
########################################
## print ###############################
########################################
########################################

print.turbo <- function(x, ...) {
	## x = object of class 'turbo'
	if(all(!is.na(x$errors))) {
		cat("All algorithms failed: \n")
		cat(error(x))
	} else {
		tab <- with(x, data.frame(method, value.objfn, itr, fpeval, objfeval, convergence, elapsed.time=runtime[,"elapsed"], row.names=seq_along(x$method)))
		tab <- tab[!x$fail,]
		print(tab, ...)
		if(any(!is.na(x$errors))) cat("\nAcceleration scheme", paste(paste(seq_along(x$method)[x$fail], " (", x$method[x$fail], ")", sep=""), collapse=", "), "failed\n")
	}
}

########################################
########################################
## pars ################################
########################################
########################################
	
pars <- function(x, ...) {
	UseMethod("pars")
}
pars.turbo <- function(x, ...) {
	## x = object of class 'turbo'
	if(length(x$method)==1) {
		return(as.vector(x$par))
	} else {
		ret <- x$pars
		rownames(ret) <- x$method
		colnames(ret) <- paste("p", 1:ncol(ret), sep="")
		return(ret[!x$fail,])
		if(any(!is.na(x$errors))) cat("\nAcceleration scheme", paste(seq_along(x$method)[x$fail], collapse=", "), "failed\n")
	}
}

########################################
########################################
## error ###############################
########################################
########################################
	
error <- function(x, ...) {
	UseMethod("error")
}
error.turbo <- function(x, ...) {
	nmethod <- length(x$method)
	which.error <- which(!is.na(x$errors))
	if(length(which.error) == 0) {
		print("There were no errors")
	} else {
		for(i in which.error) {
			##cat(paste("method ", i, " (", x$method[i], "):", sep=""), "\n", x$errors[i], "\n")
			cat(paste("method ", i, " (", x$method[i], "):", sep=""), x$errors[i])
		}
	}
}

########################################
########################################
## plot ################################
########################################
########################################

plot.turbo <- function(x, which.methods = seq_along(x$method), method.names = x$method[which.methods], xlim, ylim, ...) {
	## x = object of class 'turbo'
	## which.methods = vector identifying the subset of algorithms whose results will be plotted
	## method.names = names of the algorithms identified by 'which.methods'
	if(!x$control.run$keep.objfval)
		stop("plot methods only defined when control.run$keep.objfval=TRUE")
	trace <- x$trace.objfval[which.methods]
	select <- which((!x$fail)[which.methods])

	times <- lapply(select, function(u) c(trace[[u]]$time.before.iter["elapsed"], trace[[u]]$time.before.iter["elapsed"] + 1:x$itr[u]*trace[[u]]$time.per.iter["elapsed"]))

	max.time <- max(sapply(select, function(u) max(times[[u]])))
	lower <- max(sapply(select, function(u) min(times[[u]])))
	if(missing(xlim)) xlim <- c(lower, max.time)
	if(missing(ylim)) ylim <- range(unlist(sapply(select, function(u) trace[[u]]$trace[times[[u]] >= xlim[1] & times[[u]] <= xlim[2]])))
	plot(unlist(sapply(select, function(u) times[[u]])), unlist(sapply(select, function(u) trace[[u]]$trace)), type="n", xlab="run time (sec.)", ylab="Objective function value", main="Trace of Objective Function Value", xlim=xlim, ylim=ylim, ...)
	for(k in seq_along(select)) {
		lines(times[[select[k]]], trace[[select[k]]]$trace, col=k)
	}
	legend("topright", as.character(method.names[select]), col=seq_along(select), lwd=1)
}
	
########################################
########################################
## gradient ############################
########################################
########################################

##setGeneric("grad")
grad.numDeriv <- grad
grad <- function(x, ...) {
	UseMethod("grad")
}
grad.turbo <- function(x, objfn=x$objfn, which.methods = seq_along(x$method), method.names = x$method[which.methods], ...) {
	## x = object of class 'turbo'
	## objfn = objective function to be minimized
	## which.methods = vector identifying the subset of algorithms whose results will be plotted
	## method.names = names of the algorithms identified by 'which.methods'
	if(is.null(objfn))
		stop("objfn must be provided to compute gradient")
	subs <- (!x$fail)[which.methods]
	select.methods <- which.methods[subs]
	mat <- matrix(NA, length(select.methods), ncol(x$par))
	rownames(mat) <- method.names[subs]
	for(k in seq_along(select.methods)) {
		mat[k,] <- grad.numDeriv(objfn, x$par[select.methods[k],], method="Richardson", method.args=list(r=2), ...)
	}
	return(mat)
}

########################################
########################################
## hessian #############################
########################################
########################################

hessian.numDeriv <- hessian
hessian <- function(x, ...) {
	UseMethod("hessian")
}
hessian.turbo <- function(x, objfn=x$objfn, which.methods = seq_along(x$method), method.names = x$method[which.methods], ...) {
	## x = object of class 'turbo'
	## objfn = objective function to be minimized
	## which.methods = vector identifying the subset of algorithms whose results will be plotted
	## method.names = names of the algorithms identified by 'which.methods'
	if(is.null(objfn))
		stop("objfn must be provided to compute hessian")
	subs <- (!x$fail)[which.methods]
	select.methods <- which.methods[subs]
	lst <- vector("list", length(select.methods))
	names(lst) <- method.names[subs]
	for(k in seq_along(select.methods)) {
		lst[[k]] <- hessian.numDeriv(objfn, x$par[select.methods[k],], method="Richardson", method.args=list(r=2), ...)
	}
	return(lst)
}

########################################
########################################
## sderror #############################
########################################
########################################

stderror <- function(x, ...) {
	UseMethod("stderror")
}
stderror.turbo <- function(x, objfn=x$objfn, which.methods = seq_along(x$method), method.names = x$method[which.methods], ...) {
	## x = object of class 'turbo'
	## objfn = objective function to be minimized
	## which.methods = vector identifying the subset of algorithms whose results will be plotted
	## method.names = names of the algorithms identified by 'which.methods'
	if(is.null(objfn))
		stop("objfn must be provided to compute standard error")
	hesses <- hessian(x, objfn=objfn, which.methods=which.methods, method.names=method.names)
	ret <- t(sapply(seq_along(hesses), function(u) sqrt(diag(solve(hesses[[u]])))))
	rownames(ret) <- names(hesses)
	return(ret)
}