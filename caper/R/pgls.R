## RESTRUCTURE AND EXPANSION/MERGING OF PGLM CODE

pgls <- function(formula, data, lambda = 1.0, kappa = 1.0,  delta= 1.0, 
	             param.CI = 0.95, control = list(fnscale=-1), 
                 bounds = NULL) {

	## bounds go singular: bounds = list(delta = c(1e-04, 3), lambda = c(1e-04,  0.99999), kappa = c(1e-04, 3))
	
	## pgls replaces lik.lambda - exactly the same as a null model
	
	## all the internal functions that were here are now farmed out to externally accessible functions
	## - except Dfun because I don't know what it does!
	
	Dfun <- function(Cmat) {
		iCmat <- solve(Cmat,  tol = .Machine$double.eps)
		svdCmat <- La.svd(iCmat)
		D <- svdCmat$u %*% diag(sqrt( svdCmat$d )) %*% t(svdCmat$v)
		return( t(D) )
	}

	## think about this - allow old data + V use?
	if(! inherits(data, 'comparative.data')) stop("data is not a 'comparative' data object.")
	dname <- deparse(substitute(data))
	call <- match.call()
	
	## check for missing data in the formula and replace the data object with a complete version
	miss <- model.frame(formula, data$data, na.action=na.pass)
	miss.na <- apply(miss, 1, function(X) (any(is.na(X))))
	if(any(miss.na)) {
		miss.names <- data$phy$tip.label[miss.na]
		data <- data[-which(miss.na),]
	}
	
	# Get the design matrix, number of parameters and response variable	
	m <- model.frame(formula, data$data)
	y <- m[,1]
	x <- model.matrix(formula, m)
	k <- ncol(x)
	namey <- names(m)[1]

	# test for variables with no variance in model matrices (thx to Sarah Dryhurst)
	# - will cause singularity - lm() filters for this and aliases variables
	#   but here we'll just fail for the time being
	xVar <- apply(x,2,var)[-1] # drop intercept
	badCols <- xVar < .Machine$double.eps
	if(any(badCols)) stop('Model matrix contains columns with zero variance: ', paste(names(xVar)[badCols], collapse=', '))

	## if the comparative data doesn't contain a VCV,
	## then get one and put it in the data object too. Bit wasteful
	if(is.null(data$vcv)){
		V <- if(kappa == 1) {VCV.array(data$phy)} else {VCV.array(data$phy, dim=3)}
		data$vcv <- V
	} else {
		V <- data$vcv
	}

	## sort out the data
	nm <- names(data$data)
	n <- nrow(data$data)
	
	# if a ci is specified, check (early) that it is sensible for use at the end!
	# ha! first planned use of sequential 'or'  operator
	if(! is.null(param.CI)){
		if(! is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) 
			stop('param.CI is not a number between 0 and 1.')
	}
	
	# check and insert elements from user bounds
	usrBounds <- bounds
	bounds <- list(kappa=c(1e-6,3), lambda=c(1e-6,1), delta=c(1e-6,3))
	if(! is.null(usrBounds)){
		if(! is.list(usrBounds))
			stop("Bounds must be a list of named bounds for any or all of kappa, lambda and delta")

		usrNames <- names(usrBounds)
		badNames <- setdiff(usrNames, c('kappa', 'lambda', 'delta')) 
		if(length(badNames) > 0)
			stop("The list of bounds contains names other than kappa, lambda and delta")
		
		for(nm in usrNames){bounds[nm] <- usrBounds[nm]}
	}
	
	## check the branch length transformations to be applied
	## - gather into a named list: names are used throughout to 
	##   get the right values in the right place.
	parVals <- list(kappa=kappa, lambda=lambda, delta=delta)
	
	## - test the bounds and parameter values are sensible
	for(i in seq_along(parVals)){

		## is the parameter a single number or 'ML'
		p <- parVals[[i]]
		nm <- names(parVals)[i]

		if(length(p) > 1) stop(nm, " not of length one.")
		if(is.character(p) & p != "ML") stop(nm, " is character and not 'ML'.")

		## are the bounds of length 2, numeric and positive or zero
		bnds <- bounds[[nm]]
		if(length(bnds) > 2) stop("Bounds specified for ",nm, " not of length one.")
		if(! is.numeric(bnds)) stop("Non-numeric bounds specified for ",nm, ".")
		if(any(bnds < 0)) stop("Negative values in bounds specified for ",nm, ".")
		lb <- bnds[1]
		ub <- bnds[2]
		if(lb > ub) stop("Lower bound greater than upper bound for ",nm, ".")
		
		## are specified transforms (not 'ML') in range (also filters out negative transforms) 
		if(is.numeric(p) & ( p < lb | p > ub))
			stop(sprintf("%s value (%0.2f) is out of specified bounds [%0.2f, %0.2f]", nm, p, lb, ub))
	}
	
	if(kappa != 1 && length(dim(V)) != 3) stop("3D VCV.array needed for kappa transformation.")

	## which are being optimised
	mlVals <- sapply(parVals,  "==", "ML")

	## if any are being optimised then run pgls.likelihood as a simple optimising function,
	## returning the logLik for a particular set of transformations
	##  - start the search for ML estimates from the midpoint of the specified bounds
	
	if(any(mlVals)){
	    
    	# isolate parameters to be optimized and set to a sensible start.
    	parVals[mlVals] <- lapply(bounds, mean)[mlVals]
		# collapse list down to a vector
    	parVals <- as.numeric(parVals)
    	names(parVals) <- c("kappa", "lambda","delta")
		
		# split them up
    	optimPar <- parVals[mlVals]
    	fixedPar <- parVals[!mlVals]
    	
    	# define the optimization bounds
    	lower.b <- sapply(bounds,  "[", 1)[mlVals]
    	upper.b <- sapply(bounds,  "[", 2)[mlVals]
    	
		## TODO - could isolate single optimisations here to use optimise() rather than optim()
		## likelihood function swapped out for externally visible one
    	optim.param.vals <- optim(optimPar, fn = pgls.likelihood, # function and start vals
    	    method="L-BFGS-B", control=control, upper=upper.b, lower=lower.b, # optim control
    	    V = V, y=y, x=x, fixedPar = fixedPar, optim.output=TRUE) # arguments to function
	    
    	if(optim.param.vals$convergence != "0"){
    		stop("Problem with optim:", optim.param.vals$convergence, 
    		    optim.param.vals$message)}
	    
    	fixedPar <- c(optim.param.vals$par, fixedPar)
    	fixedPar <- fixedPar[c("kappa","lambda","delta")]
    } else {
		## reduce the list of parameters to a vector
        fixedPar <- as.numeric(parVals)
		names(fixedPar) <- c("kappa", "lambda","delta")
    }
    
	## run the likelihood function again with the fixed parameter values
	## ll <- log.likelihood(optimPar=NULL, fixedPar=fixedPar, y, x, V, optim=FALSE)
	ll <- pgls.likelihood(optimPar=NULL, fixedPar=fixedPar, y, x, V, optim.output=FALSE)
	
	## store the log likelihood of the optimized solution for use in ci.searchs
	log.lik <- ll$ll
	
	## get the transformed vcv matrix for the fitted model for use
	## in calculating the remaining outputs.
	Vt <- pgls.blenTransform(V, fixedPar)
	
	## start collating outputs:
	
	## AIC
	aic <- -2 * log.lik + 2 * k
	aicc <- -2 * log.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
	
	## coefficients
	coeffs <- ll$mu
	names(coeffs) <- colnames(x)
	varNames <- names(m)

	## predicted values
	pred <- x %*% ll$mu 
	
	##residuals
	res <- y - pred
	D <- Dfun(Vt)
	pres <- D %*% res # TODO - what is this exactly
	
	## fitted model
	fm <- list(coef = coeffs, aic = aic, log.lik = log.lik)
	
	## various variances
	RMS <- ll$s2
	RSSQ <- ll$s2 * (n - k)
	
	## null model
	xdummy <- matrix(rep(1, length(y)))
	nullMod <- pgls.likelihood(optimPar=NULL, fixedPar=fixedPar, y, xdummy, V, optim.output=FALSE)
	NMS <- nullMod$s2
	NSSQ <- nullMod$s2 * (n -1) 
	
	# Bits for parameter errors	
	errMat <- t(x)%*% solve(Vt) %*% x  
	errMat <- solve(errMat) * RMS[1] 
	sterr <- diag(errMat)
	sterr <- sqrt(sterr)
	
	
	RET <- list(model = fm, formula = formula, call=call, RMS = RMS, NMS = NMS,
	            NSSQ = NSSQ[1], RSSQ = RSSQ[1], aic = aic, aicc = aicc, n = n, k = k,
	            sterr = sterr, fitted = pred, residuals = res, phyres = pres,
	            x = x, data = data,  varNames = varNames, y = y, param = fixedPar, mlVals=mlVals,
	            namey = namey, bounds=bounds, Vt=Vt, dname=dname)
	
	class(RET) <- "pgls"
	
	## missing data
	if(any(miss.na)){
		RET$na.action <- structure(which(miss.na), class='omit', .Names=miss.names)
		
	}
	## if requested, get the confidence intervals on the optimized parameters
	## if any are actually optimised
	if(! is.null(param.CI) && any(mlVals)){
		
		## Loop over optimized parameters
		param.CI.list <- list(kappa=NULL, lambda=NULL, delta=NULL)
		mlNames <- names(mlVals)[which(mlVals)]
		
		for(param in mlNames){
			param.CI.list[[param]] <- pgls.confint(RET, param, param.CI)
		}
		
		RET$param.CI <- param.CI.list
	}
	
	return(RET)
	
}

pgls.profile <- function(pgls, which=c('lambda','kappa','delta'), N=50, param.CI=NULL){
	
	## takes a pgls model and profiles one of the branch length transformations
	
	# get the x sequence for the parameter
	which <- match.arg(which)
	bnds <- pgls$bounds[[which]]
	x <- seq(from=bnds[1], to=bnds[2], length=N)
	
	# get a matrix of parameter values
	pars <- matrix(pgls$param, nrow=N, ncol=3, byrow=TRUE)
	colnames(pars) <- names(pgls$param)
	pars[,which] <- x
	
	## now get the sequence of likelihoods for the parameter in question
	logLik <- sapply(seq_along(x), function(X){ pgls.likelihood(optimPar=NULL, fixedPar=pars[X,], y=pgls$y, x=pgls$x, V=pgls$data$vcv, optim.output=TRUE)})
	
	RET <- list(x=x,logLik=logLik, which=which, pars=pgls$param, dname=pgls$dname, formula=pgls$formula)
	class(RET) <- 'pgls.profile'
	
	# test for existing parameter ci otherwise create if asked
	if(! is.null(pgls$param.CI[which])){
		RET$ci <- pgls$param.CI[[which]]
	} else if(! is.null(param.CI)){
		RET$ci <- pgls.confint(pgls, which, param.CI)
	} 

	return(RET)
}

plot.pgls.profile <- function(x, ...){
	
	xlab <- as.expression(x$which)
	xsub <- sprintf('Data: %s; Model: %s\nkappa %0.2f; lambda %0.2f; delta %0.2f', 
	                x$dname, deparse(x$formula), x$pars['kappa'], x$pars['lambda'], x$pars['delta'])
	
	with(x, plot(logLik ~ x, type='l', xlab=xlab, ...))
	title(sub=xsub, cex.sub=0.7, line=par('mgp')[1]+1.5)
	
	
	if(! is.null(x$ci)){
		abline(v=x$ci$opt, col='red', ...)
		abline(v=x$ci$ci.val, lty=2, col='red', ...)
	}

}

pgls.confint <- function(pgls, which=c('lambda','kappa','delta'), param.CI=0.95){
	
	# Are we dealing with a same confidence interval
	# ha! first planned use of sequential 'or'  operator
	if(! is.numeric(param.CI) || param.CI <= 0 || param.CI > 1) 
		stop('ci is not a number between 0 and 1.')
	
	# find the parameter being checked
	which <- match.arg(which)
	
	# is the value in the object for this parameter an ML value?
	# - if not, then this needs to be estimated in order 
	#   to get confidence intervals.
	# - currently, bail out but could refit to model to get this.
	if(pgls$mlVals[which] == FALSE) stop('The pgls object contains a fixed, not ML, estimate of ', which)
	ML <- pgls$model$log.lik
	
	# separate out the values held constant and varied
	fix <- pgls$param
	whichNum <- which(names(fix) == which)
	opt <- fix[whichNum]
	fix <- fix[-whichNum]
	
	# only one optimPar so get bounds and two intervals
	bounds  <- pgls$bounds[[which]]
	belowML <- c(bounds[1], opt)
	aboveML <- c(opt, bounds[2])
	
	# the offset needed to get the root of the ML surface
	# at zero is  - (observed ML) + a chisq component
	
	MLdelta <- (qchisq(param.CI, 1)/2)
	offset <- (- ML) + MLdelta

	## get the model components
	y <- pgls$y
	x <- pgls$x
	V <- pgls$data$vcv
	
	## find the confidence intervals on the parameter
	## - first need to find the logLik at the bounds
	## - as long as the bound is outside the CI, can then use uniroot 
	##   to find the actual confidence interval.

	lowerBound.ll <- pgls.likelihood(structure(bounds[1], names=which), fix, y, x, V, optim.output=TRUE)
	upperBound.ll <- pgls.likelihood(structure(bounds[2], names=which), fix, y, x, V, optim.output=TRUE)
	
	lrt0 <- 2 * (ML - lowerBound.ll)
	lrt1 <- 2 * (ML - upperBound.ll)
	lowerBound.p <- 1 - pchisq(lrt0, 1)
	upperBound.p <- 1 - pchisq(lrt1, 1)
	
	## - a problem with uniroot is that the identity of the variables gets stripped
	##   which is why pgls.likelihood now has an optim.names option used here.
	ll.fun <- function(opt){
        pg <- pgls.likelihood(opt, fix, y, x, V, optim.output=TRUE, names.optim=which)
        ll <- pg + offset
        return(ll)
    }
	
	lowerCI <- if(lowerBound.ll < (ML -MLdelta)) uniroot(ll.fun , interval=belowML)$root else NA
	upperCI <- if(upperBound.ll < (ML -MLdelta)) uniroot(ll.fun , interval=aboveML)$root else NA

	return(list(opt=opt, bounds.val=bounds, bounds.p=c(lowerBound.p, upperBound.p), ci.val=c(lowerCI, upperCI), ci=param.CI))

}

pgls.likelihood <- function(optimPar, fixedPar, y, x, V, optim.output=TRUE, names.optim=NULL) {
    
	# Full ML estimation for given x and V: modified to also act as an engine for optim
	# - this is why the branch length  parameters are passed as two chunks, so that
	#   the first acts as the targets for optimisation.
	# - the function is passed named vectors containing kappa, lambda and delta
	#   which might be available for optimization (optimPar) or user defined (fixedPar)

    # merge the values of KLD from the two parameter vectors
	# if names.optim is provided then add it (uniroot in the ci.search strips it out)
	
	# Estimates the GLS parameters for given data
	get.coeffs <- function(Y, iV, X) {
		xVix <- crossprod(X, iV %*% X)
		xViy <- crossprod(X, iV %*% Y)
		mu <- solve(xVix, tol = .Machine$double.eps) %*% xViy 	#This is  a bad thing to do!!!!
		return(mu)
	}

	# Estimates the variance of a given trait (accounting for phylogeny)
	est.var <- function(y, iV, x, mu ) {
		e <- y - x %*% mu
		s2 <- crossprod(e, iV %*% e)
		n <- length(y) 
		k <- length(x[1,])
		return( s2 / (n- k) )
	}

	if(! is.null(names.optim)) names(optimPar) <- names.optim
    allPar <- c(optimPar, fixedPar)
    
	# get the transformed VCV matrix and its inverse
    V <- pgls.blenTransform(V, allPar)
	iV <- solve(V, tol = .Machine$double.eps)
	
	mu <- get.coeffs(y, iV, x)
	s2 <- est.var(y, iV, x, mu)
	n <- nrow(x)
	k <- ncol(x)
	logDetV <- determinant(V, logarithm = TRUE)$modulus[1]
	
	## Likelihood calculation
	## original pgls3.2: ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - 1)/2.0
	## Richard Duncan's implementation: log.lkhood.y <- log((1/((2*pi*sigma.sqr)^(ntax/2))) * det(VCV)^-0.5
	##  * exp(-(1 / (2 * sigma.sqr)) * t(response - design.matrix %*% alpha)
	##  %*% solve(VCV) %*% (response - design.matrix %*% alpha)))
	## ll <- log((1/((2*pi*s2)^(n/2))) * det(V) ^-0.5 
	##       * exp(-(1 / (2 * s2)) * t(y - x %*% mu)
	##       %*% solve(V) %*% (y - x %*% mu)))
	ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log( (n - k) * s2 / n) - logDetV / 2.0  - n / 2.0
	
	# if being used for optimization, only return the log likelihood
	if(optim.output) return(ll)  else return( list(ll = ll, mu = mu, s2 = s2) )
}

pgls.blenTransform <- function(V, fixedPar){
	## applies the three branch length transformations to a VCV matrix
	
    # apply transformations
	if(! is.null(fixedPar["kappa"]) && fixedPar["kappa"] != 1){
		if(length(dim(V)) < 3){
			stop('Kappa transformation requires a 3 dimensional VCV array.')
		}
	}
	
    if(fixedPar["kappa"] == 0) V <- (V > 0) else V <-  V ^ fixedPar["kappa"] # kappa catching NA^0=1
    V <- apply(V, c(1,2), sum, na.rm=TRUE) # collapse 3D array
    V <- ifelse(upper.tri(V)+lower.tri(V), V * fixedPar["lambda"], V) # lambda
    if(fixedPar["delta"] == 0) V <- (V > 0) else V <-  V ^ fixedPar["delta"] # delta catching NA^0=1

	attr(V, 'blenTransform') <- fixedPar
	return(V)
}

plot.pgls <- function(x, ...) {
	
	# layout(matrix(c(1,2,3,4), 2, 2, byrow = FALSE))
	res <- residuals(x, phylo = TRUE)
	res <- res / sqrt(var(res))[1]
	# truehist(res, xlab = "Residual value (corrected for phylogeny)")
	plot(density(res))
	qqnorm(res)
	qqline(res)
	plot(fitted(x), res, xlab = "Fitted value", ylab = "Residual value (corrected for phylogeny)"  )
	plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value")
}

summary.pgls <- function(object,...) {
	
	## call and return object
	ans <- list(call=object$call)
	class(ans) <- 'summary.pgls'
	
	## model size
	p <- object$k
	n <- object$n
	rdf <- n - p
	ans$df <- c(p, rdf)
	
	## residuals and residual standard error
	r <- object$phyres
	rss <- object$RSSQ
	resvar <- rss/rdf
	ans$sigma <- sqrt(resvar)
	ans$residuals <- r
	
	## coefficient matrix
	cf <- object$model$coef
	se <- object$sterr
	t  <- cf/se
	
	coef <- cbind(cf,se,t, 2 * ( 1 - pt(abs(t), rdf)))
	colnames(coef) <- c('Estimate','Std. Error','t value','Pr(>|t|)')
	ans$coefficients <- coef
	
	## parameter matrix
	ans$param <- object$param
	ans$mlVals <- object$mlVals
	
	if(! is.null(object$param.CI)) ans$param.CI <- object$param.CI
	
	if(! is.null(object$na.action)) ans$na.action <- object$na.action
	
	## model statistics: p includes the intercept - it is the number of columns of the design matrix
	ans$fstatistic <- c(value= ((object$NSSQ - object$RSSQ) / object$RMS) / (object$k - 1),  numdf= p-1,dendf=rdf) 
	ans$r.squared <- (object$NSSQ - object$RSSQ) / object$NSSQ
    ans$adj.r.squared <- (object$NMS - object$RMS) / object$NMS
   	
	return(ans)
	
}

print.summary.pgls <- function(x, digits = max(3, getOption("digits") - 3), ...){
	
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")	

	r <- zapsmall(quantile(x$resid), digits + 1)
	names(r) <- c("Min", "1Q", "Median", "3Q", "Max")
	cat('Residuals:\n')
	print(r, digits=digits)


	cat("\nBranch length transformations:\n\n")
	for(p in names(x$param)){
		cat(sprintf("%-6s [%s]  : %0.3f\n", p, ifelse(x$mlVals[p], " ML", "Fix"), x$param[p]))
		if(! is.null(x$param.CI[[p]])){
			blopt <- x$param.CI[[p]]
			cat(sprintf("   lower bound : %0.3f, p = %-5s\n", blopt$bounds.val[1], format.pval(blopt$bounds.p[1])))
			cat(sprintf("   upper bound : %0.3f, p = %-5s\n", blopt$bounds.val[2], format.pval(blopt$bounds.p[2])))
			cat(sprintf("   %2.1f%% CI   : (%0.3f, %0.3f)\n", blopt$ci *100,blopt$ci.val[1], blopt$ci.val[2]))
		}
	}


	cat('\nCoefficients:\n')
	printCoefmat(x$coef)

    cat("\nResidual standard error:", format(signif(x$sigma, 
        digits)), "on", x$df[2L], "degrees of freedom\n")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
	cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
    cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, 
        digits = digits), "\nF-statistic:", formatC(x$fstatistic[1L], 
        digits = digits), "on", x$fstatistic[2L], "and", 
        x$fstatistic[3L], "DF,  p-value:", format.pval(pf(x$fstatistic[1L], 
            x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE), 
            digits = digits), "\n")
    

}

print.pgls <- function(x,  digits = max(3, getOption("digits") - 3), ...){
	
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
        quote = FALSE)
	cat("\n")
}

coef.pgls <- function(object, ...){
	
	cf <- object$model$coef
	nm <- rownames(cf)
	cf <- structure(as.vector(cf), names=nm)
	return(cf)
	
}

# This returns the residuals from the model
## CDLO - argument name changed for consistency with S3 generic
residuals.pgls <- function(object, phylo = FALSE, ...) {
    ret <- NULL
	if(phylo == FALSE){ret <- object$res} else {ret <- object$phyres}
	return(ret)
}

# This returns the fitted values
## CDLO - argument name changed for consistency with S3 generic
fitted.pgls <- function(object, ...){
    ret <- object$fitted
    return(ret)
}

# This predicts for given x
## CDLO - argument name changed for consistency with S3 generic
## CDLO - argument name of x changed to discriminate from generic to plot and print

## predict.pgls <- function(object, pred.x, ...) {
##     mu <- as.matrix(coef(object) )
##     ret <- cbind(1,  pred.x) %*% t(mu)
##     return(ret)
## }

## - rewritten by CDLO to use standard method as in lm - prompted by Mike Steiper
predict.pgls <- function (object, newdata=NULL, ...){
	
    # pull the data from the model if no new data is provided
    if(is.null(newdata)){
    	newdata <- object$data$data
    }
    
    # turn that into a design matrix
    # need to drop the response from the formula
    dmat <- model.matrix(delete.response(terms(formula(object))), data=newdata)
    
    # multiply through by the coefficients
    mu <- as.matrix(coef(object))
    ret <- dmat %*% mu
    return(ret)
}


## enables the generic AIC methods for objects and lists of objects 
logLik.pgls <- function(object, REML = FALSE, ...){
	
	val <- object$model$log.lik
	
	attr(val, "nall") <- object$n
    attr(val, "nobs") <- object$n
    attr(val, "df") <- object$k
    class(val) <- "logLik"
    val
}

## generic function for number of observations. 
nobs.pgls <- function(object, ...) length(resid(object))


## # This returns the AICc
## ## CDLO - argument name changed for consistency with S3 generic
## AICc.pgls <- function(object) {
##     ret <- object$aicc
##     return(ret[1])
## }

anova.pgls <- function(object, ...){
	
	## SEQUENTIAL SUMS OF SQUARES.
	## ASSUMES ORDER OF TERMS PRESERVE MARGINALITY
	
    if (length(list(object, ...)) > 1L){
        return(anova.pglslist(object, ...))
	} else {
	    data <- object$data
	    tlabels <- attr( terms(object$formula), "term.labels")
		k <- object$k
		n <- object$n
		NR <- length(tlabels) + 1
	
		# track residual ss and residual df and get residuals and df of null model 
		rss <- resdf <- rep(NA, NR)
		rss[1] <- object$NSSQ
		resdf[1] <- n - 1
	
		lm <- object$param['lambda']
		dl <- object$param['delta']
		kp <- object$param['kappa']
	
		# fit the sequential models
	    for( i in 1:(k-1)) {
	    		fmla <- as.formula(paste( object$namey, " ~ ", paste(tlabels[1:i], collapse = "+") ))
	    		plm <- pgls(fmla, data, lambda=lm, delta=dl, kappa=kp)
	    		rss[i+1] <- plm$RSSQ
	    		resdf[i+1] <- (n - 1) - plm$k + 1
	    }

		ss <- c(abs(diff(rss)), object$RSSQ)
		df <- c(abs(diff(resdf)), n -k)
		ms <- ss/df
		fval <- ms / ms[NR]
	    P <-  pf(fval, df, df[NR], lower.tail=FALSE)

	    table <- data.frame(df, ss, ms, f=fval, P)
	    table[length(P), 4:5] <- NA
	    dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", 
	        "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
	    #if (attr(object$terms, "intercept")) 
	    #    table <- table[-1, ]
	    structure(table, heading = c("Analysis of Variance Table", 
	        sprintf("Sequential SS for pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n", lm, dl,kp), 
	        paste("Response:", deparse(formula(object)[[2L]]))), 
	        class = c("anova", "data.frame"))
	}
}

anova.pglslist <- function(object, ..., scale = 0, test = "F"){
	
    objects <- list(object, ...)

	## check the models use the same response
    responses <- as.character(lapply(objects, function(x) deparse(terms(x$formula)[[2L]])))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
        objects <- objects[sameresp]
        warning("models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from ", "model 1")
    }
    
	## check the models have the same number of cases (not actually that they are the same values)
    ns <- sapply(objects, function(x) length(x$residuals))
    if (any(ns != ns[1L])) 
        stop("models were not all fitted to the same size of dataset")
    
	## check that the model parameters are the same
	param <- sapply(objects, '[[', 'param')
	paramChk <- apply(param, 1, function(X) all(X == X[1]))
	if(! all(paramChk))
	    stop('models were fitted with different branch length transformations.')
	
    nmodels <- length(objects)
    if (nmodels == 1) 
        return(anova(object))
    resdf <- as.numeric(lapply(objects, function(X) X$n - X$k))
    resdev <- as.numeric(lapply(objects, '[[', 'RSSQ'))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
        -diff(resdev)))
    variables <- lapply(objects, function(x) paste(deparse(formula(x)), 
        collapse = "\n"))
    dimnames(table) <- list(1L:nmodels, c("Res.Df", "RSS", "Df", 
        "Sum of Sq"))
    title <- "Analysis of Variance Table"
    subtitle <- sprintf("pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n", 
                        param['lambda', 1], param['delta', 1], param['kappa', 1])
    topnote <- paste("Model ", format(1L:nmodels), ": ", variables, 
        sep = "", collapse = "\n")
    if (!is.null(test)) {
        bigmodel <- order(resdf)[1L]
        scale <- if (scale > 0) 
            scale
        else resdev[bigmodel]/resdf[bigmodel]
        table <- stat.anova(table = table, test = test, scale = scale, 
            df.scale = resdf[bigmodel], n = length(objects[bigmodel$residuals]))
    }
    structure(table, heading = c(title, subtitle, topnote), class = c("anova", 
        "data.frame"))
}

