cov_sp_arg_check <- function(coords, sp.type, sp.par, 
	error.var, smoothness, finescale.var, pcoords, 
	D, Dp, Dop)
{
	#check coords argument
	if(!is.numeric(coords) || !is.matrix(coords)){ stop("coords must be a numeric matrix") }

	#check sp.type arguments
	if(!valid_sp_type(sp.type)){ stop("sp.type is not a valid covariance type") }

	#check sp.par argument
	if(!is.numeric(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(!is.vector(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(length(sp.par) != 2){ stop("sp.par must be a numeric vector of length 2") }
	if(!(min(sp.par) > 0)){ stop("sp.par must have positive elements") }

	#check error.var argument
	if(!is.numeric(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(!is.vector(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(length(error.var) > 1){ stop("error.var must be a numeric vector of length 1") }
	if(min(error.var) < 0){ stop("error.var must be non-negative") }

	#check smoothness argument
	if(!is.numeric(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(!is.vector(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(length(smoothness) > 1){ stop("smoothness must be a numeric vector of length 1") }
	if(min(smoothness) <= 0){ stop("smoothness must be positive") }

	#check finescale.var argument
	if(!is.numeric(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(!is.vector(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(length(finescale.var) > 1){ stop("finescale.var must be a numeric vector of length 1") }
	if(min(finescale.var) < 0){ stop("finescale.var must be non-negative") }

	#check pcoords argument
	if(!is.null(pcoords))
	{
		if(!is.numeric(pcoords) || !is.matrix(pcoords)){ stop("pcoords must be a numeric matrix") }
		if(ncol(coords) != ncol(pcoords))
		{ 
			stop("coords and pcoords must have the same number of columns") 
		}
	}
	
	#check D argument
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D) ){ stop("D must be a numeric matrix if provided") } 
		if(min(D) < 0){ stop("D cannot have negative elements") }
		if(!isSymmetric(D)){ stop("D must be a symmetric matrix") } 
		if(nrow(D) != nrow(coords))
		{
			stop("D must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Dp argument
	if(!is.null(Dp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dp is provided") }
		if(!is.numeric(Dp) || !is.matrix(Dp)){ stop("Dp must be a numeric matrix if provided") } 
		if(min(Dp) < 0){ stop("Dp cannot have negative elements") }
		if(!isSymmetric(Dp)){ stop("Dp must be a symmetric matrix") }
		if(nrow(Dp) != nrow(pcoords))
		{
			stop("Dp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Dp argument
	if(!is.null(Dop))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dop is provided") }
		if(!is.numeric(Dop)  || !is.matrix(Dop)){ stop("Dop must be a numeric matrix if provided") } 
		if(min(Dop) < 0){ stop("Dop cannot have negative elements") }
		if(ncol(Dop) != nrow(pcoords))
		{
			stop("Dop must have ncols equal to nrows of pcoords")
		}
		if(nrow(Dop) != nrow(coords))
		{
			stop("Dop must have nrows equal to nrows of coords")
		}
	}
}

valid_sp_type <- function(sp.type)
{
	#returns TRUE if sp.type is equal to the below options
	#otherwise it returns FALSE
	return((sp.type == "exponential" ||
		sp.type == "gaussian" ||
		sp.type == "matern" ||
		sp.type == "spherical")||
		sp.type == "matern2")
}

valid_t_type <- function(t.type)
{
	#returns TRUE if sp.type is equal to the below options
	#otherwise it returns FALSE
	return((t.type == "ar1"))
}

valid_decomp_type <- function(method)
{
	return((method == "eigen" || method == "chol" || method == "svd"))
}

cov_st_arg_check <- function(coords, time, sp.type, sp.par, 
	error.var, smoothness, finescale.var, t.type, t.par, pcoords, 
	ptime, D, Dp, Dop, T, Tp, Top)
{
	#check coords argument
	if(!is.numeric(coords) || !is.matrix(coords)){ stop("coords must be a numeric matrix") }

	#check time argument
	if(!is.numeric(time) || nrow(time) != nrow(coords))
	{ 
		stop("time must be a numeric matrix with nrows equal to nrow(coords) (or a vector with length equal to nrow(coords))")
	} 

	#check sp.type arguments
	if(!valid_sp_type(sp.type)){ stop("specified sp.type is not a valid covariance type") }

	#check sp.par argument
	if(!is.numeric(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(!is.vector(sp.par)){ stop("sp.par must be a numeric vector of length 2") }
	if(length(sp.par) != 2){ stop("sp.par must be a numeric vector of length 2") }
	if(!(min(sp.par) > 0)){ stop("sp.par must have positive elements") }


	#check error.var argument
	if(!is.numeric(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(!is.vector(error.var)){ stop("error.var must be a numeric vector of length 1") }
	if(length(error.var) > 1){ stop("error.var must be a numeric vector of length 1") }
	if(min(error.var) < 0){ stop("error.var must be non-negative") }


	#check smoothness argument
	if(!is.numeric(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(!is.vector(smoothness)){ stop("smoothness must be a numeric vector of length 1") }
	if(length(smoothness) > 1){ stop("smoothness must be a numeric vector of length 1") }
	if(min(smoothness) <= 0){ stop("smoothness must be positive") }

	#check finescale.var argument
	if(!is.numeric(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(!is.vector(finescale.var)){ stop("finescale.var must be a numeric vector of length 1") }
	if(length(finescale.var) > 1){ stop("finescale.var must be a numeric vector of length 1") }
	if(min(finescale.var) < 0){ stop("finescale.var must be non-negative") }

	#check t.type argument
	if(!valid_t_type(t.type)){ stop("specified t.type does not match available options") }

	#check t.par argument
	if(!is.numeric(t.par) || t.par < 0 || t.par >= 1){ stop("t.par must be in range [0, 1)") }
	
	#check pcoords argument
	if(!is.null(pcoords))
	{
		if(!is.numeric(pcoords)){ stop("pcoords must be a numeric matrix") }
		if(ncol(coords) != ncol(pcoords))
		{ 
			stop("coords and pcoords must have the same number of columns") 
		}
		if(is.null(ptime)){ stop("ptime must be supplied when pcoords is supplied.") }
		
		#check ptime argument
		if(!is.numeric(ptime) || nrow(ptime) != nrow(pcoords))
		{ 
			stop("ptime must be a numeric matrix with nrows equal to nrow(pcoords) (or a vector with length equal to nrow(pcoords)")
		} 
	}

	
	#check D argument
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D) ){ stop("D must be a numeric matrix if provided") } 
		if(min(D) < 0){ stop("D cannot have negative elements") }
		if(!isSymmetric(D)){ stop("D must be a symmetric matrix") } 
		if(nrow(D) != nrow(coords))
		{
			stop("D must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Dp argument
	if(!is.null(Dp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dp is provided") }
		if(!is.numeric(Dp) || !is.matrix(Dp)){ stop("Dp must be a numeric matrix if provided") } 
		if(min(Dp) < 0){ stop("Dp cannot have negative elements") }
		if(!isSymmetric(Dp)){ stop("Dp must be a symmetric matrix") }
		if(nrow(Dp) != nrow(pcoords))
		{
			stop("Dp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Dp argument
	if(!is.null(Dop))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Dop is provided") }
		if(!is.numeric(Dop)  || !is.matrix(Dop)){ stop("Dop must be a numeric matrix if provided") } 
		if(min(Dop) < 0){ stop("Dop cannot have negative elements") }
		if(ncol(Dop) != nrow(pcoords))
		{
			stop("Dop must have ncols equal to nrows of pcoords")
		}
		if(nrow(Dop) != nrow(coords))
		{
			stop("Dop must have nrows equal to nrows of coords")
		}
	}

	#check T argument
	if(!is.null(T))
	{
		if(!is.numeric(T) || !is.matrix(T) ){ stop("T must be a numeric matrix if provided") } 
		if(min(T) < 0){ stop("T cannot have negative elements") }
		if(!isSymmetric(T)){ stop("T must be a symmetric matrix") } 
		if(nrow(T) != nrow(coords))
		{
			stop("T must be a matrix with nrows and ncols equal to nrows of coords")
		}
	}

	#check Tp argument
	if(!is.null(Tp))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Tp is provided") }
		if(!is.numeric(Tp) || !is.matrix(Tp)){ stop("Tp must be a numeric matrix if provided") } 
		if(min(Tp) < 0){ stop("Tp cannot have negative elements") }
		if(!isSymmetric(Tp)){ stop("Tp must be a symmetric matrix") }
		if(nrow(Tp) != nrow(pcoords))
		{
			stop("Tp must be a matrix with nrows and ncols equal to nrows of pcoords")
		}
	}

	#check Top argument
	if(!is.null(Top))
	{
		if(is.null(pcoords)){ stop("pcoords must be supplied when Top is provided") }
		if(!is.numeric(Top)  || !is.matrix(Top)){ stop("Top must be a numeric matrix if provided") } 
		if(min(Top) < 0){ stop("Top cannot have negative elements") }
		if(ncol(Top) != nrow(pcoords))
		{
			stop("Top must have ncols equal to nrows of pcoords")
		}
		if(nrow(Top) != nrow(coords))
		{
			stop("Top must have nrows equal to nrows of coords")
		}
	}
}

decomp_cov_check_arg <- function(V, method, checkSymmetric = TRUE)
{
	if(!is.matrix(V) || !is.numeric(V))
	{
		stop("V must be a numeric matrix")
	}
	#Removed because sometimes a symmetric matrix may not be due to numerical imprecision
	#if(checkSymmetric)
	#{
	#	if(!isSymmetric(V)){ stop("V must be symmetric") }
	#}
	if(!(method == "eigen" || method == "chol" || method == "svd"))
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}
}

maxlik_cov_sp_check_arg <- function(X, y, coords, sp.type, 
	range.par, error.ratio, smoothness, D, reml, lower, upper)
{
	if(!is.numeric(X) || !is.matrix(X)){ stop("X must be a numeric matrix") }
	if(!is.numeric(y)){ stop("y must be numeric") }
	if(!is.numeric(coords) || !is.numeric(coords)){ stop("coords must be a numeric matrix") }
	if(!valid_sp_type(sp.type)){ stop("specified sp.type is not a valid covariance type") }
	if(!(range.par > 0)){ stop("range.par must be positive") }
	if(!(error.ratio >= 0)){ stop("range.par must be non-negative") }
	if(!(smoothness > 0)){ stop("smoothness must be positive") }
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D))
		{ 
			stop("If supplied, D must be a numeric matrix") 
		}
	}
	if(!is.logical(reml)){ stop("reml must be a logical value") }
	if(!is.null(lower))
	{
		if(!is.numeric(lower)){ stop("lower must be a numeric vector")} 
	}
	if(!is.null(upper))
	{
		if(!is.numeric(upper)){ stop("upper must be a numeric vector")} 
	}
	if(!is.null(lower) && !is.null(upper))
	{
		if(length(lower) != length(upper)){ stop("lower and upper should have the same length") }
	}
}

maxlik_cov_st_check_arg <- function(X, y, coords, time = time, sp.type, 
	range.par, error.ratio, smoothness, t.type = "ar1", D, T, reml, lower, upper)
{
	if(!is.numeric(X) || !is.matrix(X)){ stop("X must be a numeric matrix") }
	if(!is.numeric(y)){ stop("y must be numeric") }
	if(!is.numeric(coords) || !is.numeric(coords)){ stop("coords must be a numeric matrix") }
	if(!valid_sp_type(sp.type)){ stop("specified sp.type is not a valid covariance type") }
	if(!(range.par > 0)){ stop("range.par must be positive") }
	if(!(error.ratio >= 0)){ stop("range.par must be non-negative") }
	if(!(smoothness > 0)){ stop("smoothness must be positive") }
	if(t.type != "ar1"){ stop("ar1 is the only valid option for t.type") }
	if(!is.null(D))
	{
		if(!is.numeric(D) || !is.matrix(D))
		{ 
			stop("If supplied, D must be a numeric matrix") 
		}
	}
	if(!is.null(T))
	{
		if(!is.numeric(T) || !is.matrix(T))
		{ 
			stop("If supplied, T must be a numeric matrix") 
		}
	}

	if(!is.logical(reml)){ stop("reml must be a logical value") }
	if(!is.null(lower))
	{
		if(!is.numeric(lower)){ stop("lower must be a numeric vector")} 
	}
	if(!is.null(upper))
	{
		if(!is.numeric(upper)){ stop("upper must be a numeric vector")} 
	}
	if(!is.null(lower) && !is.null(upper))
	{
		if(length(lower) != length(upper)){ stop("lower and upper should have the same length") }
	}
}

krige_arg_check <- function(y, V, Vp, Vop, X, Xp, m, 
	nsim, Ve.diag, method)
{

	if(!is.numeric(y))
	{
		stop("y must be numeric")
	}
	if(!is.vector(y))
	{
		stop("y must be a vector")
	}

	n <- length(y)

	if(!is.matrix(V) || !is.numeric(V) || (nrow(V)!= ncol(V)))
	{
		stop("V must be a square numeric matrix")
	}
	if(nrow(V) != n)
	{
		stop("nrow(V) must equal length(y))")
	}

	if(!is.matrix(Vp) || !is.numeric(Vp) || (nrow(Vp)!= ncol(Vp)))
	{
		stop("Vp must be a square numeric matrix")
	}

	if(!is.matrix(Vop) || !is.numeric(Vop))
	{
		stop("Vop must be a numeric matrix")
	}

	if(n != nrow(Vop))
	{
		stop("length(y) must equal nrow(Vop)")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncol(Vp) must equal ncol(Vop)")
	}
	
	if((is.null(X) && !is.null(Xp)) || (!is.null(X) && is.null(Xp)))
	{
		stop("If X is supplied, Xp must also be supplied (and vice versa)")
	}

	if(!is.null(X))
	{
		if(nrow(X) != n)
		{
			stop("nrow(X) must equal length(y)")
		}
		if(nrow(Xp) != nrow(Vp))
		{
			stop("nrow(Xp) must equal nrow(Vp)")
		}
		if(ncol(Xp) != ncol(X))
		{
			stop("ncol(Xp) must equal ncol(X)")
		}
	}
	
	if(!(is.numeric(m) && length(m) == 1))
	{
		stop("m must be a numeric vector of length 1")
	}

	if(!(is.numeric(nsim) && length(nsim) == 1 && (nsim >= 0)))
	{
		stop("nsim must be a nonnegative value")
	}

	if(nsim >= 1)
	{
		if(!(is.vector(Ve.diag) && length(Ve.diag) == n && min(Ve.diag) >= 0))
		{
			stop("Ve.diag must be a vector of length(y) with non-negative values")
		}
		if(!(method == "eigen" || method == "chol" || method == "svd"))
		{
			stop("method must be 'eigen', 'chol', or 'svd'")
		}
	}
	else
	{
		Ve.diag <- 1
	}
	
	# change method to number for C++ function
	if(method == "eigen"){ method <- 1 }
	else if(method == "chol"){ method <- 2}
	else{ method <- 3 }
	
	return(list(method = method, Ve.diag = Ve.diag))
}

condnorm_par_arg_check <- function(y, V, Vp, Vop, coeff, X, Xp, method)
{
	n <- length(y)
	nk <- length(coeff)

	if(!is.numeric(y))
	{
		stop("y must be a numeric vector")
	}
	if(!is.matrix(V) || !is.numeric(V) || (nrow(V)!= ncol(V)))
	{
		stop("V must be a square numeric matrix")
	}
	if(!is.matrix(Vp) || !is.numeric(Vp) || (nrow(Vp)!= ncol(Vp)))
	{
		stop("Vp must be a square numeric matrix")
	}
	if(!is.matrix(Vop) || !is.numeric(Vop))
	{
		stop("Vop must be a numeric matrix")
	}
	if(length(y) != nrow(V))
	{
		stop("length of y must match nrows of V")
	}
	if(length(y) != nrow(Vop))
	{
		stop("length of y must match nrows of Vop")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncols of Vp must match ncols of Vop")
	}

	if(!is.numeric(coeff))
	{
		stop("coeff must be a numeric vector")
	}
	if((is.null(X) && !is.null(Xp)) || (!is.null(X) && is.null(Xp)))
	{
		stop("If X is supplied, Xp must also be supplied (and vice versa)")
	}
	if(!is.null(X))
	{
		if(nrow(X) != n)
		{
			stop("nrows of X must match length of y")
		}
		if(ncol(X) != nk)
		{
			stop("ncols of X must match length of coeff")
		}
		if(nrow(Xp) != nrow(Vp))
		{
			stop("nrows of Xp must match nrows of Vp")
		}
		if(ncol(Xp) != ncol(X))
		{
			stop("ncols of Xp must match ncols of X")
		}
	}
	if(!valid_decomp_type)
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}
}

rmvnorm_arg_check <- function(nsim, mu, V, method)
{
	if(!is.numeric(nsim) || !is.numeric(mu) || ! is.numeric(V))
	{
		stop("nsim, mu, and V arguments must all be numeric")
	}
	#if(!isSymmetric(V) || !is.matrix(V))
	if(!is.matrix(V))
	{
		stop("V must be a symmetric matrix")
	}
	if(length(mu) != nrow(V))
	{
		stop("The length of mu must equal nrows of V")
	}
	if(!(method == "eigen" || method == "chol" || method == "svd"))
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}
}

rcondsim_arg_check <- function(nsim, y, V, Vp, Vop, Ve.diag, method = "eigen", krige.obj)
{
	if(!(is.numeric(nsim)))
	{
		stop("nsim must be a positive number")
	}
	if(length(nsim) != 1)
	{
		stop("nsim must be a vector of length 1")
	}
	if(!(nsim > 1))
	{
		stop("nsim must be at least 1")
	}
	if(!is.numeric(y) || !is.vector(y))
	{
		stop("y must be a numeric vector")
	}
	if(!is.matrix(V) || !is.numeric(V) || (nrow(V)!= ncol(V)))
	{
		stop("V must be a square numeric matrix")
	}
	if(length(y) != nrow(V))
	{
		stop("length(y) must equal nrow(V)")
	}
	if(!is.matrix(Vp) || !is.numeric(Vp) || (nrow(Vp)!= ncol(Vp)))
	{
		stop("Vp must be a square numeric matrix")
	}
	if(!is.matrix(Vop) || !is.numeric(Vop))
	{
		stop("Vop must be a numeric matrix")
	}
	if(length(y) != nrow(Vop))
	{
		stop("length(y) must equal nrow(Vop)")
	}
	if(ncol(Vp) != ncol(Vop))
	{
		stop("ncol(Vp) must equal ncol(Vop)")
	}
	if(!(is.vector(Ve.diag) && (is.numeric(Ve.diag))))
	{
		stop("Ve.diag must be a numeric vector")
	}
	if(length(y) != length(Ve.diag))
	{
		stop("Ve.diag must have the same length as y")
	}
	if(!(method == "eigen" || method == "chol" || method == "svd"))
	{
		stop("method must be 'eigen', 'chol', or 'svd'")
	}	
	if(!is.list(krige.obj))
	{
		stop("krige.obj should be an object returned by a kriging function and should contain w, a matrix of prediction weights for the observed data.")		
	}
	if(is.null(krige.obj$w))
	{
		stop("krige.obj should be an object returned by a kriging function and should contain w, a matrix of prediction weights for the observed data.")		
	}
}
