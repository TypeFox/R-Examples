#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

.rho <- function(x, lamb, derive = 0, type = c("EL", "ET", "CUE"), k = 1)
	{

	type <- match.arg(type)
	lamb <- matrix(lamb, ncol = 1)
	gml <- x%*%lamb*k
	if (derive == 0)
		{
		if (type == "EL")
			{
			if (any(gml>=1))
				stop("Computation of Lambda fails because NAs produced by log(1-gt*l)")
			rhomat <- log(1 - gml) 
			}
		if (type == "ET")
			rhomat <- -exp(gml)
		if (type == "CUE")
			rhomat <- -gml -0.5*gml^2
		}
	if (derive==1)
		{
		if (type == "EL")
			rhomat <- -1/(1 - gml) 
		if (type == "ET")
			rhomat <- -exp(gml)
		if (type == "CUE")
			rhomat <- -1 - gml
		}
	if (derive==2)
		{
		if (type == "EL")
			rhomat <- -1/(1 - gml)^2 
			
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -rep(1,nrow(x))
		}
	return(c(rhomat))
	}

.getCgelLam <- function(gt, l0, type = c('EL', 'ET', 'CUE'), method = c("nlminb", "optim", "constrOptim"), control=list(), k = 1, alpha = 0.01)
	{
	type <- match.arg(type)
	method <- match.arg(method)
	fct <- function(l, X)
		{
		r1 <- colMeans(.rho(gt,l,derive=1,type=type,k=k)*X)
		crossprod(r1) + alpha*crossprod(l)
		}
	Dfct <- function(l, X)
		{
		r2 <- .rho(X,l,derive=2,type=type,k=k)
		r1 <- .rho(X,l,derive=1,type=type,k=k)
		H <- t(X*r2)%*%X/nrow(X)
		2*H%*%colMeans(r1*X) + 2*alpha*l
		}

	if (method == "nlminb")
		res <- nlminb(l0, fct, Dfct, X = gt, control = control) 
	if (method == "optim")
		res <- optim(l0, fct, Dfct, X = gt, method="BFGS", control = control) 
	if (method == "constrOptim")
		{
		ci <- -rep(1,nrow(gt))
		res <- constrOptim(rep(0,ncol(gt)),fct,Dfct,-gt,ci,control=control,X=gt)
		}
	if (method == "optim")
		{
		conv <- list(convergence = res$convergence, counts = res$counts, message = res$message)
	} else {
		conv <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
		}

	return(list(lambda = res$par, convergence = conv, 
			obj = mean(.rho(gt,res$par, derive=0,type=type,k=k))))	
	}

getLamb <- function(gt, l0, type = c('EL', 'ET', 'CUE', "ETEL"), tol_lam = 1e-7, maxiterlam = 100, tol_obj = 1e-7, 
		k = 1, 	method = c("nlminb", "optim", "iter"), control=list())
	{
	method <- match.arg(method)
	if(is.null(dim(gt)))
		gt <- matrix(gt,ncol=1)

	if (method == "iter")
		{
		if (type == "ETEL")
			type = "ET"
		for (i in 1:maxiterlam)
			{
			r1 <- .rho(gt,l0,derive=1,type=type,k=k)
			r2 <- .rho(gt,l0,derive=2,type=type,k=k)
			F <- -colMeans(r1*gt)
			J <- crossprod(r2*gt,gt)
			if (sum(abs(F))<tol_obj)
				{
				conv <- list(convergence="Tolerance for the FOC reached")
				break
				}
			P <- solve(J,F)
			if (sum(abs(P))<tol_lam)
				{
				conv <- list(convergence="Tolerance on lambda reached")	
				break
				}
			l0 <- l0 + P
			conv <- list(convergence="maxiterlam reached")
			}
		}
	 else
		{
		fct <- function(l,X)
			{
			r0 <- .rho(X,l,derive=0,type=type,k=k)
			-mean(r0)
			}
		Dfct <- function(l,X)
			{
			r1 <- .rho(X,l,derive=1,type=type,k=k)
		        -colMeans(r1*X)
			}
		DDfct <- function(l,X)
			{
			r2 <- .rho(X,l,derive=2,type=type,k=k)
			-t(X*r2)%*%X/nrow(X)
			}
		if (type == "ETEL")
			{
			type = "ET"
			ci <- -rep(1,nrow(gt))
			res <- constrOptim(rep(0,ncol(gt)),fct,Dfct,-gt,ci,control=control,X=gt)
			}
		else
			{
			if (method=="optim")
				{
				if (type != "EL")
					res <- optim(rep(0,ncol(gt)),fct,gr=Dfct,X=gt,method="BFGS",control=control)
				else
					{		
					ci <- -rep(1,nrow(gt))
					res <- constrOptim(rep(0,ncol(gt)),fct,Dfct,-gt,ci,control=control,X=gt)
					}
				}
			else
				res <- nlminb(rep(0,ncol(gt)), fct, gradient = Dfct, hessian = DDfct, X = gt, control = control)
			}

		l0 <- res$par
		if (method == "optim" | method == "constrOptim")
			conv <- list(convergence = res$convergence, counts = res$counts, message = res$message)
		if(method == "nlminb")
			conv <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
		}
	return(list(lambda = l0, convergence = conv, obj = mean(.rho(gt,l0,derive=0,type=type,k=k))))
	}

smoothG <- function (x, bw = bwAndrews, prewhite = 1, ar.method = "ols", weights = weightsAndrews,
			kernel = c("Bartlett", "Parzen", "Truncated", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
			tol = 1e-7) 
	{
	kernel <- match.arg(kernel)
	approx <- match.arg(approx)
		
	n <- nrow(x)
	if (is.function(weights))
		{
                        class(x) <- "gmmFct"
			w <- weights(x, bw = bw, kernel = kernel,  
			prewhite = prewhite, ar.method = ar.method, tol = tol, 
			verbose = FALSE, approx = approx)
		}
		else
			w <- weights


	rt <- length(w)
	if (rt >= 2)
		{
		rt <- length(w)
		if (rt>1)
			{
			w <- c(w[rt:2], w)
			w <- w / sum(w)
			w <- kernel(w[rt:length(w)])
			}
		else
			w <- kernel(1)

		x <- kernapply(x,w)		
		sx <- list("smoothx" = x, "kern_weights" = w)
		return(sx)		
		}
	else
		sx <- list("smoothx" = x,"kern_weights" = kernel(1))
		return(sx)		
	}


gel <- function(g, x, tet0, gradv = NULL, smooth = FALSE, type = c("EL", "ET", "CUE", "ETEL"), 
                kernel = c("Truncated", "Bartlett"), bw = bwAndrews, approx = c("AR(1)", 
    		"ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol_weights = 1e-7, tol_lam = 1e-9, tol_obj = 1e-9, 
		tol_mom = 1e-9, maxiterlam = 100, constraint = FALSE, optfct = c("optim", "optimize", "nlminb"), 
                optlam = c("nlminb", "optim", "iter"), data, Lambdacontrol = list(), model = TRUE, X = FALSE, Y = FALSE, TypeGel = "baseGel", alpha = NULL, ...)
	{

	type <- match.arg(type)
	optfct <- match.arg(optfct)
	optlam <- match.arg(optlam)
	weights <- weightsAndrews
	approx <- match.arg(approx)
	kernel <- match.arg(kernel)
	if(missing(data))
		data<-NULL
	all_args <- list(g = g, x = x, tet0 = tet0, gradv = gradv, smooth = smooth, type = type,
                kernel = kernel, bw = bw, approx = approx, prewhite = prewhite, ar.method = ar.method, 
		tol_weights = tol_weights, tol_lam = tol_lam, tol_obj = tol_obj, tol_mom = tol_mom, 
		maxiterlam = maxiterlam, constraint = constraint, optfct = optfct, weights = weights,
                optlam = optlam, model = model, X = X, Y = Y, TypeGel = TypeGel, call = match.call(), 
		Lambdacontrol = Lambdacontrol, alpha = alpha, data = data)

	class(all_args)<-TypeGel
	Model_info<-getModel(all_args)
	z <- momentEstim(Model_info, ...)

	class(z) <- "gel"
	return(z)
	}


.thetf <- function(tet, P, output=c("obj","all"), l0Env)
    {
    output <- match.arg(output)
    gt <- P$g(tet, P$dat)
    l0 <- get("l0",envir=l0Env)
    if (is.null(P$CGEL))
	{
	if (P$optlam != "optim" & P$type == "EL") 
	    {
	    lamb <- try(getLamb(gt, l0, type = P$type, tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, 
			tol_obj = P$tol_obj, k = P$k1/P$k2, control = P$Lambdacontrol, 
			method = P$optlam), silent = TRUE)
            if(class(lamb) == "try-error")
		    lamb <- getLamb(gt, l0, type = P$type, tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, 
			tol_obj = P$tol_obj, k = P$k1/P$k2, control = P$Lambdacontrol, method = "optim")
	    }
	    else
		    lamb <- getLamb(gt, l0, type = P$type,	 tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, 
			tol_obj = P$tol_obj, k = P$k1/P$k2, control = P$Lambdacontrol, method = P$optlam)
	}
    else
	{
	if (P$type=="ETEL")
		stop("CGEL not implemented for ETEL")
	lamb <- try(.getCgelLam(gt, l0, type = P$type, method = "nlminb", control=P$Lambdacontrol, 
			k = P$k1/P$k2, alpha = P$CGEL),silent=TRUE)
	if (class(lamb) == "try-error")
		lamb <- try(.getCgelLam(gt, l0, type = P$type, method = "constrOptim", control=P$Lambdacontrol, 
			k = P$k1/P$k2, alpha = P$CGEL),silent=TRUE)
	}
	
    obj <- mean(.rho(gt, lamb$lambda, type = P$typet, derive = 0, k = P$k1/P$k2))
    assign("l0",lamb$lambda,envir=l0Env)
    if(output == "obj")
	    return(obj)
    else
	    return(list(obj = obj, lambda = lamb, gt = gt))
    }


