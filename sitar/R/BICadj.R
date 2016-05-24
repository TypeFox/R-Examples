	BICadj <- function(..., pattern=NULL)
#	input SITAR model(s)
#	returns BIC for power transformed y
#	looks for log, sqrt or ^ - assumes no multipliers
{	ARG <- match.call(expand.dots=FALSE)$...
	if (!is.null(pattern)) pattern <- ls(envir=parent.frame(), pattern=pattern)
	ARG <- unique(c(unlist(sapply(ARG, deparse)), pattern))
	dev <- sapply(ARG, function(obj) {
		obj <- get(obj)
		if (is.character(try(ll <- logLik(obj), TRUE))) return(NA)
#	check for call.sitar$y or else call$y or else call$formula
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
#	check for call$y or call$model
		if (!is.null(obj$call$y)) ycall <- obj$call$y 
			else if (!is.null(obj$call$model)) ycall <- obj$call$model[[2]] 
			else if (!is.null(obj$call$formula)) ycall <- obj$call$formula[[2]] 
			else return(NA)
		lambda <- 1
		if (length(ycall) == 1) yt <- ycall else {
			yt <- ycall[[2]]
			fun <- ycall[1]
			if (fun == "log()") lambda <- 0 else
			if (fun == "sqrt()") lambda <- 0.5 else 
			if (fun == "`^`()") lambda <- eval(ycall[[3]])
		}
		y <- eval(yt, eval(obj$call$data, sys.parent()))
		sly <- ifelse(lambda == 1, 0, sum(log(y)))
		BIC(ll) - 2 * ((lambda - 1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))		
	})
	dev <- dev[!is.na(dev)]
	if (length(dev) == 0) return(invisible())
	round(dev[order(dev)], 1)
}

	AICadj <- function(..., k=2, pattern=NULL)
#	input SITAR model(s)
#	returns AIC for power transformed y
#	looks for log, sqrt or ^ - assumes no multipliers
{	ARG <- match.call(expand.dots=FALSE)$...
	if (!is.null(pattern)) pattern <- ls(envir=parent.frame(), pattern=pattern)
	ARG <- unique(c(unlist(sapply(ARG, deparse)), pattern))
	dev <- sapply(ARG, function(obj) {
		obj <- get(obj)
		if (is.character(try(ll <- logLik(obj), TRUE))) return(NA)
#	check for call.sitar$y or else call$y or else call$formula
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
#	check for call$y or call$model
		if (!is.null(obj$call$y)) ycall <- obj$call$y 
			else if (!is.null(obj$call$model)) ycall <- obj$call$model[[2]] 
			else if (!is.null(obj$call$formula)) ycall <- obj$call$formula[[2]] 
			else return(NA)
		lambda <- 1
		if (length(ycall) == 1) yt <- ycall else {
			yt <- ycall[[2]]
			fun <- ycall[1]
			if (fun == "log()") lambda <- 0 else
			if (fun == "sqrt()") lambda <- 0.5 else 
			if (fun == "`^`()") lambda <- eval(ycall[[3]])
		}
		y <- eval(yt, eval(obj$call$data, sys.parent()))
		sly <- ifelse(lambda == 1, 0, sum(log(y)))
		AIC(ll, k=k) - 2 * ((lambda - 1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))		
	})
	dev <- dev[!is.na(dev)]
	if (length(dev) == 0) return(invisible())
	round(dev[order(dev)], 1)
}
