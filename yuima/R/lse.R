##estimate theta2 by LSE
##function name LSE1


lse <- function(yuima, start, lower, upper, method="BFGS", ...){ 	

	call <- match.call()
	
	if( missing(yuima))
	yuima.stop("yuima object is missing.")
	
## param handling
	
## FIXME: maybe we should choose initial values at random within lower/upper
##        at present, qmle stops	
	if( missing(start) ) 
	yuima.stop("Starting values for the parameters are missing.")
	
	diff.par <- yuima@model@parameter@diffusion
	drift.par <- yuima@model@parameter@drift
	jump.par <- yuima@model@parameter@jump
	measure.par <- yuima@model@parameter@measure
	common.par <- yuima@model@parameter@common

	fullcoef <- c(diff.par, drift.par)
	npar <- length(fullcoef)
	
		
	nm <- names(start)
    oo <- match(nm, fullcoef)
    if(any(is.na(oo))) 
	yuima.stop("some named arguments in 'start' are not arguments to the supplied yuima model")
    start <- start[order(oo)]
    nm <- names(start)
	
	
	
	idx.diff <- match(diff.par, nm)
	idx.drift <- match(drift.par, nm)

	
	tmplower <- as.list( rep( -Inf, length(nm)))
	names(tmplower) <- nm	
	if(!missing(lower)){
		idx <- match(names(lower), names(tmplower))
		if(any(is.na(idx)))
		yuima.stop("names in 'lower' do not match names fo parameters")
		tmplower[ idx ] <- lower	
	}
	lower <- tmplower
	
	tmpupper <- as.list( rep( Inf, length(nm)))
	names(tmpupper) <- nm	
	if(!missing(upper)){
		idx <- match(names(upper), names(tmpupper))
		if(any(is.na(idx)))
		yuima.stop("names in 'lower' do not match names fo parameters")
		tmpupper[ idx ] <- upper	
	}
	upper <- tmpupper

	
	
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	for(t in 1:(n-1))
	 env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	
##objective function
	f <-function(theta){
		names(theta) <- drift.par
		tmp <- env$deltaX - env$h * drift.term(yuima, theta, env)[-n,]
		ret <- t(tmp) %*% tmp
		return(sum(ret))
	}

	mydots <- as.list(call)[-(1:2)]
	mydots$fixed <- NULL
	mydots$fn <- as.name("f")
	mydots$start <- NULL
	mydots$par <- unlist(start)
	mydots$hessian <- FALSE
	mydots$upper <- unlist( upper[ nm[idx.diff] ])
	mydots$lower <- unlist( lower[ nm[idx.diff] ])
	
	
	if(length(start)>1){ #multidimensional optim				
		oout <- do.call(optim, args=mydots)
	} else { ### one dimensional optim
		mydots$f <- mydots$fn
		mydots$fn <- NULL
		mydots$par <- NULL
		mydots$hessian <- NULL	
		mydots$method <- NULL	
		mydots$interval <- as.numeric(c(lower[drift.par],upper[drift.par])) 
		mydots$lower <- NULL	
		mydots$upper <- NULL
		opt1 <- do.call(optimize, args=mydots)
#opt1 <- optimize(f, ...) ## an interval should be provided
		oout <- list(par = opt1$minimum, value = opt1$objective)
	} 
	
	return(oout$par)            
}





