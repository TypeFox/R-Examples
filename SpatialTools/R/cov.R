cov.sp <- function(coords, sp.type = "exponential", 
	sp.par = stop("specify sp.par argument"), 
	error.var = 0, smoothness = 0.5, finescale.var = 0, pcoords = NULL, 
	D = NULL, Dp = NULL, Dop = NULL)
{
	#check validity of arguments
	cov_sp_arg_check(coords, sp.type, sp.par, error.var, smoothness, 
		finescale.var, pcoords, D, Dp, Dop)

	# calculate covariance matrix for observed responses
	# calculate distance matrix if not supplied
	if(is.null(D))
	{ 
		V <- simple.cov.sp(dist1(coords), sp.type, sp.par, error.var, smoothness, finescale.var)
	}
	else
	{
		V <- simple.cov.sp(D, sp.type, sp.par, error.var, smoothness, finescale.var)
	}
	
	#calculate covariance matrices for predicted responses and
	#between observed responses and predicted responses
	#calculate Dp and Dop if necessary
	if(!is.null(pcoords))
	{
		if(is.null(Dp))
		{ 
			Vp <- simple.cov.sp(dist1(pcoords), sp.type, sp.par, error.var = 0, smoothness, finescale.var) 
		}
		else
		{
			Vp <- simple.cov.sp(Dp, sp.type, sp.par, error.var = 0, smoothness, finescale.var) 
		}
		if(is.null(Dop))
		{
			Vop <- simple.cov.sp(dist2(coords, pcoords), sp.type, sp.par, error.var = 0, smoothness, finescale.var)
		}
		else
		{
			Vop <- simple.cov.sp(Dop, sp.type, sp.par, error.var = 0, smoothness, finescale.var)
		}
	}
	
	if(is.null(pcoords))
	{
		return(list(V = V)) 
	}else
	{
		return(list(V = V, Vp = Vp, Vop = Vop)) 
	}
}

simple.cov.sp <- function(D, sp.type, sp.par, error.var, smoothness, finescale.var)
{
	if(sp.type == "exponential")
	{
		V <- sp.par[1]*exp(-D/sp.par[2])

	}else if(sp.type == "gaussian")
	{
		V <- sp.par[1]*exp(-(D/sp.par[2])^2)

	}else if(sp.type == "matern")
	{
		sD <- D/sp.par[2]
		V <- (D > 0) * sp.par[1]*(2^(1-smoothness)/gamma(smoothness)*sD^smoothness*besselK(sD, nu = smoothness))
		V[is.nan(V)] <- sp.par[1]
	}else if(sp.type == "matern2")
	{
		sD <- 2 * sqrt(smoothness) * D/sp.par[2]
		V <- (D > 0) * sp.par[1]*(2^(1-smoothness)/gamma(smoothness)*sD^smoothness*besselK(sD, nu = smoothness))
		V[is.nan(V)] <- sp.par[1]	
	}else if(sp.type == "spherical")
	{
		sD <- D/sp.par[2]
		V <- sp.par[1]*(1 - 1.5*sD + 0.5*(sD)^3)*(D < sp.par[2])
	}
	return(V + (D == 0)*(finescale.var + error.var))
}

cov.st <- function(coords, time, sp.type = "exponential", sp.par = stop("specify sp.par argument"), 
	error.var = 0, smoothness = 0.5, finescale.var = 0, t.type = "ar1", t.par = .5, pcoords = NULL, 
	ptime = NULL, D = NULL, Dp = NULL, Dop = NULL, T = NULL, Tp = NULL, Top = NULL)
{

	if(!is.matrix(time)){ time <- matrix(time) }
	if(!is.null(ptime))
	{
		if(!is.matrix(ptime)){ ptime <- matrix(ptime) }
	}

	cov_st_arg_check(coords, time, sp.type, sp.par, 
		error.var, smoothness, finescale.var, t.type, t.par, pcoords, 
		ptime, D, Dp, Dop, T, Tp, Top)

	if(is.null(D)){ D <- dist1(coords) }
	V.sp <- simple.cov.sp(D, sp.type, sp.par, error.var = 0, smoothness, finescale.var = 0)

	if(is.null(T)){ T <- dist1(matrix(time)) }
	V.time <- simple.cov.time(T, t.type, t.par)

	#Only add finescale or measurement error when at the same location and time	
	V <- V.sp * V.time + (D==0) * (T==0) * (finescale.var + error.var)

	if(!is.null(pcoords))
	{
		if(is.null(ptime)){ stop("ptime must be specified since pcoords specified") }
		if(is.null(Dp)){ Dp <- dist1(pcoords) }
		if(is.null(Dop)){ Dop <- dist2(coords, pcoords) }
		if(is.null(Tp)){ Tp <- dist1(matrix(ptime)) }
		if(is.null(Top)){ Top <- dist2(matrix(time), matrix(ptime)) }
		
		Vp.sp <- simple.cov.sp(Dp, sp.type, sp.par, error.var = 0, smoothness, finescale.var = 0)
		Vp.time <- simple.cov.time(Tp, t.type, t.par)
		Vp <- Vp.sp*Vp.time + (Dp == 0)*finescale.var

		Vop.sp <- simple.cov.sp(Dop, sp.type, sp.par, error.var = 0, smoothness, finescale.var = 0)
		Vop.time <- simple.cov.time(Top, t.type, t.par)
		Vop <- Vop.sp*Vop.time + (Dop == 0)*finescale.var
	}

	if(is.null(pcoords))
	{
		return(list(V = V)) 
	}else
	{
		return(list(V = V, Vp = Vp, Vop = Vop)) 
	}
}

simple.cov.time <- function(T, t.type, t.par)
{
	if(t.type == "ar1")
	{
		V.time <- t.par^T	
	}
	V.time
}

decomp.cov <- function(V, method = "eigen")
{
	decomp_cov_check_arg(V, method)
	
	if(method == "eigen")
	{
		eigenV <- eigen(V)
		return(eigenV$vectors %*% diag(sqrt(pmax(eigenV$values,0))))
	}else if(method == "chol")
	{ 
		return(t(chol(V))) 
	}else if(method == "svd")
	{
		svdV <- svd(V)
		return(tcrossprod(svdV$u %*% diag(sqrt(svdV$d)), svdV$v))
	}
}

#marginally faster than decomp.V for small matrices, but much slower for large
#decomp.cov2 <- function(V, method = "eigen")
#{
#	decomp_cov_check_arg(V, method)
#	
#	if(method == "eigen")
#	{
#		.Call( "decomp_cov", Vs = V, 
#			method = 1, PACKAGE = "SpatialTools")
#	}else if(method == "chol")
#	{ 
#		.Call( "decomp_cov", Vs = V, 
#			method = 2, PACKAGE = "SpatialTools")
#	}else if(method == "svd")
#	{
#		.Call( "decomp_cov", Vs = V, 
#			method = 3, PACKAGE = "SpatialTools")
#	}
#}

maxlik.cov.sp <- function(X, y, coords, sp.type = "exponential", 
	range.par = stop("specify range.par argument"), 
	error.ratio = stop("specify error.ratio argument"), smoothness = 0.5,  
	D = NULL, reml = TRUE, lower = NULL, upper = NULL, 
	control = list(TRACE = TRUE), optimizer="nlminb")
{

	y <- as.vector(y)
	maxlik_cov_sp_check_arg(X, y, coords, sp.type, 
		range.par, error.ratio, smoothness, D, reml, lower, upper)

	if(is.null(D)){ D <- dist1(coords) }

	#construct vector of initial parameters to pass to objective function
	parms <- c(range.par, error.ratio)
	if(sp.type == "matern"){ parms <- c(parms, smoothness) }

	#define lower and upper bounds if not supplied by user
	bounds <- bounds.cov.sp(sp.type)
	if(is.null(lower)){ lower <- bounds$lower.bound }
	if(is.null(upper)){ upper <- bounds$upper.bound	}

	if(optimizer=="nlminb")
	{
		min.out <- stats::nlminb(start=parms, 
			objective = logLik.cov.sp, X=X, y=y, D=D,
			sp.type = sp.type, reml=reml, minus2 = TRUE, lower=lower, upper=upper, 
			control = control)
			
		range.par <- min.out$par[1]
		error.ratio <- min.out$par[2]
		if(sp.type == "matern"){ smoothness <- min.out$par[3] }
	}	

	V <- simple.cov.sp(D = D, sp.type = sp.type, sp.par = c(1, range.par), 
		error.var = error.ratio, smoothness = smoothness, finescale.var = 0)

	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)
	fitted <- X %*% solve(XtViX, crossprod(ViX, y))
	r <- y - fitted
	n <- length(y)
	df.resid <- n - ncol(X)
	rtVir <- crossprod(r, solve(V,r))
	if(!reml)
	{
		sigmasq <- c(rtVir/n)
	}else
	{
		sigmasq <- c(rtVir/df.resid)
	}
	
	error.var <- sigmasq * error.ratio

	sp.par <- c(sigmasq, range.par)

	return(list(sp.type = sp.type, sp.par = sp.par, 
		error.var = error.var, smoothness = smoothness, par = min.out$par,
		objective = min.out$objective, convergence = min.out$convergence,
		message = min.out$message, iterations = min.out$iter, 
		evaluations = min.out$eval))
}

logLik.cov.sp <- function(par, X, y, D, sp.type, reml = FALSE, minus2 = FALSE)
{
	range.par <- par[1]
	error.ratio <- par[2]
	
	if(sp.type == "matern")
	{ 
		smoothness <- par[3]
	}else
	{
		smoothness <- 0.5
	}
	
	V <- simple.cov.sp(D = D, sp.type = sp.type, sp.par = c(1, range.par), 
		error.var = error.ratio, smoothness = smoothness, finescale.var = 0)

	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)
	fitted <- X %*% solve(XtViX, crossprod(ViX, y))
	r <- y - fitted
	n <- length(y)
	df.resid <- n - ncol(X)
	rtVir <- crossprod(r, solve(V,r))

	if(!reml)
	{
		sigmasq <- c(rtVir/n)
	
		lik <- determinant(V, logarithm=TRUE)$modulus + 
		n + length(r) * log(2*pi) + length(r)*log(sigmasq)
	}else
	{
		sigmasq <- c(rtVir/df.resid)
	
		lik <- determinant(V, logarithm=TRUE)$modulus + 
		determinant(XtViX, logarithm=TRUE)$modulus +
		(n - ncol(X)) + n * log(2*pi) + (n - ncol(X))*log(sigmasq)
	}
	if(!minus2)
	{
		lik <- -lik/2
	}

	return(lik)
}

maxlik.cov.st <- function(X, y, coords, time, sp.type = "exponential", 
	range.par = stop("specify range.par argument"), 
	error.ratio = stop("specify error.ratio argument"), 
	smoothness = 0.5, t.type = "ar1", t.par = .5, D = NULL, T = NULL, 
	reml = TRUE, lower = NULL, upper = NULL, control = list(TRACE = TRUE), optimizer="nlminb")
{
	y <- as.vector(y)
	maxlik_cov_st_check_arg(X = X, y = y, coords = coords, time = time, 
		sp.type = sp.type, range.par = range.par, error.ratio = error.ratio, 
		smoothness = smoothness, t.type = t.type, 
		D = D, T = T, reml = reml, lower = lower, upper = upper)

	if(is.null(D)){ D <- dist1(coords) }
	if(is.null(T)){ T <- dist1(matrix(time)) }

	#construct vector of initial parameters to pass to objective function
	parms <- c(range.par, error.ratio)
	if(sp.type == "matern"){ parms <- c(parms, smoothness) }

	parms <- c(parms, t.par)

	#define lower and upper bounds if not supplied by user
	bounds <- bounds.cov.sp(sp.type)
	if(is.null(lower)){ lower <- bounds$lower.bound }
	if(is.null(upper)){ upper <- bounds$upper.bound }
	lower <- c(lower, 0)
	upper <- c(upper, .999)

	if(optimizer=="nlminb")
	{
		min.out <- stats::nlminb(start=parms, 
			objective = logLik.cov.st, X=X, y=y, D=D,
			sp.type = sp.type, t.type = t.type, T=T, 
			reml=reml, minus2 = TRUE, 
			lower=lower, upper=upper, control = control)
			
		range.par <- min.out$par[1]
		error.ratio <- min.out$par[2]
		if(sp.type == "matern")
		{ 
			smoothness <- min.out$par[3]
			t.par <- min.out$par[4]
		}else
		{
			t.par <- min.out$par[3]
		}
	}	

	V.sp <- simple.cov.sp(D, sp.type, sp.par = c(1, range.par), 
		error.var = 0, smoothness, finescale.var = 0)

	V.time <- simple.cov.time(T, t.type, t.par)
	
	V <- V.sp * V.time + diag(nrow(V.sp)) * error.ratio

	n <- length(y)
	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)
	fitted <- X %*% solve(XtViX, crossprod(ViX, y))
	r <- y - fitted
	df.resid <- n - ncol(X)
	rtVir <- crossprod(r, solve(V,r))

	if(!reml)
	{
		sigmasq <- as.vector(rtVir/n)
	}else
	{
		sigmasq <- as.vector(rtVir/df.resid)
	}
	
	error.var <- sigmasq * error.ratio
	sp.par <- c(sigmasq, range.par)

	return(list(sp.type = sp.type, sp.par = sp.par, 
		error.var = error.var, t.par = t.par, 
		smoothness = smoothness, par = min.out$par,
		objective = min.out$objective, convergence = min.out$convergence,
		message = min.out$message, iterations = min.out$iter, 
		evaluations = min.out$eval))
}

logLik.cov.st <- function(par, X, y, D, T, sp.type, t.type, reml = FALSE, minus2 = FALSE)
{
	range.par <- par[1]
	error.ratio <- par[2]
	
	if(sp.type == "matern")
	{ 
		smoothness <- par[3]
		t.par <- par[4]
	}else
	{
		smoothness <- 0.5
		t.par <- par[3]
	}
	
	
	V.sp <- simple.cov.sp(D, sp.type, sp.par = c(1, range.par), 
		error.var = 0, smoothness, finescale.var = 0)

	V.time <- simple.cov.time(T, t.type, t.par)
	
	V <- V.sp * V.time + diag(nrow(V.sp)) * error.ratio

	ViX <- solve(V, X)
	XtViX <- crossprod(ViX, X)
	fitted <- X %*% solve(XtViX, crossprod(ViX, y))
	r <- y - fitted
	n <- length(y)
	df.resid <- n - ncol(X)
	rtVir <- crossprod(r, solve(V,r))

	if(!reml)
	{
		sigmasq <- c(rtVir/n)
	
		lik <- determinant(V, logarithm=TRUE)$modulus + 
		  n + length(r) * log(2*pi) + length(r)*log(sigmasq)
	}else
	{
		sigmasq <- c(rtVir/df.resid)
	
		lik <- determinant(V, logarithm=TRUE)$modulus + 
		  determinant(XtViX, logarithm=TRUE)$modulus +
		  (n - ncol(X)) + n * log(2*pi) + (n - ncol(X))*log(sigmasq)
	}
	if(!minus2)
	{
		lik <- -lik/2
	}

	return(lik)
}

bounds.cov.sp <- function(sp.type)
{
	lower <- c(.001, 0)
	upper <- c(Inf, 1)

	if(sp.type == "matern")
	{
		lower <- c(lower, .001)
		upper <- c(upper, 10)
	}
	
	return(list(lower.bound = lower, upper.bound = upper))
}
