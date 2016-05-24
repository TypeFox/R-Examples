# Author: cns
# fitting functions


dlmodeler.build.unknowns <- function(model)
{
	list.unknowns <- function(mat,name,fun,times=FALSE) {
		idx <- which(is.na(mat),arr.ind=TRUE)
		apply(idx, 1, function(x) list(name=name,i=x[1],j=x[2],fun=fun,times=times))
	}
	model$unknowns <- c(model$unknowns, list.unknowns(model$Ht,"Ht",function(x) exp(x)))
	model$unknowns <- c(model$unknowns, list.unknowns(model$Rt,"Rt",function(x) x))
	model$unknowns <- c(model$unknowns, list.unknowns(model$Qt,"Qt",function(x) exp(x)))
	model$unknowns <- c(model$unknowns, list.unknowns(model$Tt,"Tt",function(x) x))
	return(model)
}


dlmodeler.build.function <- function(model)
{
	if( class(model)!='dlmodeler' ) stop("model should be of class 'dlmodeler'")
	if( NROW(model$Zt)!=1 ) stop("multivariate case is not yet implemented yet") # TODO
	# TODO parametrization Qt = p*Ht
	# TODO call dlmodeler.build.unknowns
	
	build.fun <- function(par=rep(0,length(model$unknowns))) {
		built.model <- model
		nb.unknowns <- length(model$unknowns)
		if( nb.unknowns>0 ) for( pos in seq(1,nb.unknowns) ) {
			unknown <- model$unknowns[[pos]]
			if( unknown$name=="Ht" ) { built.model$Ht[unknown$i,unknown$j] <- unknown$fun(par[pos]) }
			else if( unknown$name=="Rt" ) { built.model$Rt[unknown$i,unknown$j] <- unknown$fun(par[pos]) }
			else if( unknown$name=="Qt" ) { built.model$Qt[unknown$i,unknown$j] <- unknown$fun(par[pos]) }
			else if( unknown$name=="Tt" ) { built.model$Tt[unknown$i,unknown$j] <- unknown$fun(par[pos]) }
			else stop("unknown: ", unknown$name)
		}
		built.model$unknowns <- NULL
		return(built.model)
	}
	
	return(build.fun)
}


AIC.dlmodeler.fit <-
		function(object, ..., k=2)
{
	if( class(object)!='dlmodeler.fit' ) stop("argument should be of class 'dlmodeler.fit'")
	# AIC, see page 152 of the Durbin & Koopman book
	ll <- object$logLik
	npar <- length(object$par) + sum(object$model$P0inf)
	return( -2*ll + k*npar )
}



logLik.dlmodeler.fit <-
		function(object, ...)
{
	if( class(object)!="dlmodeler.fit" ) stop("argument should be of class 'dlmodeler.fit'")
	if( is.null(object$filter) ) stop("log-likelihood is not available: call dlmodeler.fit() with filter=TRUE")
	if( is.null(object$filter$logLik) ) stop("log-likelihood is not available: call dlmodeler.fit() with filter=TRUE")
	# return log likelihood, if it has been computed
	return(object$filter$logLik)
}



logLik.dlmodeler.filtered <-
		function(object, ...)
{
	if( class(object)!="dlmodeler.filtered" ) stop("argument should be of class 'dlmodeler.filtered'")
	if( is.null(object$logLik) ) stop("log-likelihood is not available: call dlmodeler.filter() with logLik=TRUE")
	# return log likelihood, if it has been computed
	return(object$logLik)
}



dlmodeler.fit <-
		function(yt, model=NULL, method=c("MLE","MSE","MAD","MAPE","sigma"), ...)
{
	if( is.null(model) ) {
		if( method[1]=="MLE" ) return(dlmodeler.fit.MLE(yt, ...))
		if( method[1]=="MSE" ) return(dlmodeler.fit.MSE(yt, ...))
		if( method[1]=="MAD" ) return(dlmodeler.fit.MAD(yt, ...))
		if( method[1]=="MAPE" ) return(dlmodeler.fit.MAPE(yt, ...))
		if( method[1]=="sigma" ) return(dlmodeler.fit.sigma(yt, ...))
		stop( "method unknown: ", method[1] )
	} else {
		model <- dlmodeler.build.unknowns(model)
		if( method[1]=="MLE" ) return(dlmodeler.fit.MLE(yt, dlmodeler.build.function(model), rep(0,length(model$unknowns)), ...))
		if( method[1]=="MSE" ) return(dlmodeler.fit.MSE(yt, dlmodeler.build.function(model), rep(0,length(model$unknowns)), ...))
		if( method[1]=="MAPE" ) return(dlmodeler.fit.MAPE(yt, dlmodeler.build.function(model), rep(0,length(model$unknowns)), ...))
		if( method[1]=="MAD" ) return(dlmodeler.fit.MAD(yt, dlmodeler.build.function(model), rep(0,length(model$unknowns)), ...))
		if( method[1]=="sigma" ) return(dlmodeler.fit.sigma(yt, dlmodeler.build.function(model), rep(0,length(model$unknowns)), ...))
		stop( "method unknown: ", method[1] )
	}
}



dlmodeler.fit.MLE <-
		function(yt, build.fun, par, backend=c('KFAS','FKF','dlm'),
				method="L-BFGS-B", verbose=FALSE, silent=FALSE,
				filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.filter(yt, build.fun(p, ...), backend=backend, logLik=TRUE, filter=FALSE)
		if( verbose ) cat(p,':',-f$logLik,'\n')
		return(-f$logLik)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( !silent ) cat(opt$message,":",opt$convergence,"in",opt$counts[["function"]],"iterations\n")
	opt$model <- build.fun(opt$par, ...)
	if(filter|smooth) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- -opt$value
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}



dlmodeler.fit.MSE <-
		function(yt, build.fun, par, ahead, iters=NCOL(yt)-ahead-start-1, step=1, start=1,
				backend=c('KFAS','FKF','dlm'), method="L-BFGS-B", verbose=FALSE, silent=FALSE,
				filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.forecast(yt, build.fun(p, ...),
				ahead=ahead, iters=iters, step=step, start=start, backend=backend)
		mse <- mean((f$yhat-f$y)^2,na.rm=TRUE)
		if( verbose ) cat(p,':',mse,'\n')
		return(mse)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( !silent ) cat(opt$message,":",opt$convergence,"in",opt$counts[["function"]],"iterations\n")
	opt$model <- build.fun(opt$par, ...)
	if(filter|smooth) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}	



dlmodeler.fit.MAD <-
		function(yt, build.fun, par, ahead, iters=NCOL(yt)-ahead-start-1, step=1, start=1,
				backend=c('KFAS','FKF','dlm'), method="L-BFGS-B", verbose=FALSE, silent=FALSE,
				filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.forecast(yt, build.fun(p, ...),
				ahead=ahead, iters=iters, step=step, start=start, backend=backend)
		mse <- mean(abs(f$yhat-f$y),na.rm=TRUE)
		if( verbose ) cat(p,':',mse,'\n')
		return(mse)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( !silent ) cat(opt$message,":",opt$convergence,"in",opt$counts[["function"]],"iterations\n")
	opt$model <- build.fun(opt$par, ...)
	if(filter|smooth) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}	



dlmodeler.fit.MAPE <-
		function(yt, build.fun, par, ahead, iters=NCOL(yt)-ahead-start-1, step=1, start=1,
				backend=c('KFAS','FKF','dlm'), method="L-BFGS-B", verbose=FALSE, silent=FALSE,
				filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	if( min(yt,na.rm=TRUE)<=0 ) stop("MAPE is only valid for positive time-series")
	fit.fun <- function(p, ...) {
		f <- dlmodeler.forecast(yt, build.fun(p, ...),
				ahead=ahead, iters=iters, step=step, start=start, backend=backend)
		mse <- mean(abs(f$yhat-f$y)/f$y,na.rm=TRUE)
		if( verbose ) cat(p,':',mse,'\n')
		return(mse)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( !silent ) cat(opt$message,":",opt$convergence,"in",opt$counts[["function"]],"iterations\n")
	opt$model <- build.fun(opt$par, ...)
	if(filter|smooth) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}	



dlmodeler.fit.sigma <-
		function(yt, build.fun, par, backend=c('KFAS','FKF','dlm'), method="L-BFGS-B",
				verbose=FALSE, silent=FALSE, filter=TRUE, smooth=FALSE, raw.result=FALSE, ...)
{
	fit.fun <- function(p, ...) {
		f <- dlmodeler.filter(yt, build.fun(p, ...), backend=backend)
		sigma <- sqrt(var(f$f[1,1:NCOL(yt)] - yt[1,]))
		if( verbose ) cat(p,':',sigma,'\n')
		return(sigma)
	}
	opt <- optim(par, fit.fun, ..., method=method)
	if( !silent ) cat(opt$message,":",opt$convergence,"in",opt$counts[["function"]],"iterations\n")
	opt$model <- build.fun(opt$par, ...)
	if(filter|smooth) opt$filtered <- dlmodeler.filter(yt, opt$model, backend=backend, smooth=smooth, raw.result=raw.result, logLik=TRUE)
	opt$logLik <- opt$filtered$logLik
	opt$par0 <- par
	class(opt) <- "dlmodeler.fit"
	return(opt)
}
