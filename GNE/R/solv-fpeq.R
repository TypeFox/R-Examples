#############################################################################
#   Copyright (c) 2012 Christophe Dutang                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### fixed-point equations in GNE
###
###         R functions
### 

fpeq <- function(xinit,	fn, merit, 
	method=c("pure", "UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"), 
	control=list(), stepfunc, argstep, silent=TRUE, order.method=1, ...)
{
	method <- match.arg(method, c("pure","UR", "vH", "RRE", "MPE", "SqRRE", "SqMPE"))
	if(method %in% c("SqRRE", "SqMPE"))
	{
		order.method <- ifelse(order.method == 1, 2, order.method)
		method <- substr(method, 3, 5)
	}
	
	
	
	if(method == "UR")
	{
		if(missing(fn) || missing(xinit) || missing(stepfunc) || missing(merit))
			stop("missing parameters for UR method")
		if(!missing(argstep))
			finalstep <- function(x) stepfunc(x, argstep)
		else
			finalstep <- stepfunc
	}
	
	if(method == "vH")
		if(missing(fn) || missing(xinit) || missing(merit))
			stop("missing parameters for vH method")
	
	noitercount <- FALSE
		
	#inner iterations to compute fixed point function fn
	inner.counts.fpfn <- c(0, 0) #call to gap and grad gap
	inner.iter.fpfn <- 0 #iter number
		
	#fixed point function
	wrapfn <- function(x) 
	{
		fx <- fn(x)
		inner.counts.fpfn <<- inner.counts.fpfn + fx$counts
		inner.iter.fpfn <<- inner.iter.fpfn + fx$iter
		fx$par
	}
	
	if(!is.null(merit))
	{
		#inner iterations to compute merit function merit
		inner.counts.merit <- c(0, 0) #call to gap and grad gap
		inner.iter.merit <- 0 #iter number
		wrapmerit <- function(x, ...) 
		{
			vx <- merit(x, ...)
			inner.counts.merit <<- inner.counts.merit + vx$counts
			inner.iter.merit <<- inner.iter.merit + vx$iter
			vx$value
		}
	}else
	{
		wrapmerit <- NULL
		inner.counts.merit <- inner.iter.merit <- 0
	}	
	
	#default control parameters
	#for relaxationAlgoVH()
	con <- list(sigma=9/10, beta=1/2, tol=1e-6, maxit=100, echo=0) 
	namc <- names(con)
	con[namc <- names(control)] <- control
	#for fpiter()
	confpiter <- list(tol=1e-6, maxiter=100, trace=FALSE) 
	if(!is.null(control$tol))
		confpiter$tol <- control$tol
	if(!is.null(control$maxit))
		confpiter$maxiter <- control$maxit
	if(!is.null(control$trace))
		confpiter$trace <- control$trace
	#for squarem()
	consquarem <- list(tol=1e-6, maxiter=100, trace=FALSE, K=order.method, 
		method=ifelse(order.method==1, 1*(method == "RRE")+2*(method == "MPE"), method))
	if(!is.null(control$tol))
		consquarem$tol <- control$tol
	if(!is.null(control$maxit))
		consquarem$maxiter <- control$maxit
	if(!is.null(control$trace))
		consquarem$trace <- control$trace
	
	
	if(method == "pure" && !is.null(wrapmerit))
	{
		myGNE <- try( pureFP(xinit, wrapfn, wrapmerit, control=con, ...), 
					 silent=silent)
	}
	if(method == "pure" && is.null(wrapmerit))
	{
		myGNE <- try( fpiter(xinit, wrapfn, wrapmerit, control=confpiter, ...), 
					 silent=silent)
		
		if(class(myGNE) != "try-error")
		{
			if(!silent)
				print(myGNE)
			myGNE$value <- max(abs(fn(myGNE$par)$par))
			myGNE$counts <- c(fn=myGNE$fpevals, merit=myGNE$objfevals)
			myGNE$code <- 1*(myGNE$convergence == 0) + 4*(myGNE$fpevals[1] >= confpiter$maxiter) 
			myGNE$code <- myGNE$code + 2*(myGNE$convergence != 0 & myGNE$fpevals[1] <= consquarem$maxiter)
			
		}	
	}
	
	
	if(method == "UR")
	{
		myGNE <- try( relaxationAlgoUR(xinit, finalstep, wrapfn, wrapmerit, 
					control=con, ...), silent=silent)
	}
	
	if(method == "vH")
	{
		myGNE <- try( relaxationAlgoVH(xinit, wrapfn, wrapmerit, 
					control=con, ...), silent=silent)
	}
	
	if(method %in% c("RRE","MPE","SqRRE","SqMPE") && !is.null(wrapmerit))
	{
		myGNE <- try( extrapolFP(xinit, wrapfn, wrapmerit, 
					control=con, method=method, ...), silent=FALSE)
	}
	
	if(method %in% c("RRE","MPE","SqRRE","SqMPE") && is.null(wrapmerit))
	{
		myGNE <- try( squarem(xinit, wrapfn, control=consquarem, ...), silent=FALSE)
		
		if(class(myGNE) != "try-error")
		{
			if(!silent)
				print(myGNE)
			
			myGNE$value <- max(abs(fn(myGNE$par)$par))
			myGNE$counts <- c(fn=myGNE$fpevals, merit=myGNE$objfevals)
			myGNE$code <- 1*(myGNE$convergence == 0) + 4*(myGNE$fpevals[1] > consquarem$maxiter) 
			myGNE$code <- myGNE$code + 2*(myGNE$convergence != 0 & myGNE$fpevals[1] <= consquarem$maxiter)
		}			
	}
	
	
	if(class(myGNE) != "try-error")
		res <- list(par=myGNE$par, value=myGNE$value,
			outer.counts=myGNE$counts, outer.iter=myGNE$counts[1], 
			code=myGNE$code, inner.iter=inner.iter.fpfn+inner.iter.merit, 
			inner.counts=inner.counts.fpfn+inner.counts.merit,
			message=myGNE$message)
	else 
		res <- list(par= NA, value=NA, outer.counts=NA, outer.iter=NA, code=100, 
			message=paste("Error in the fixed-point problem:", myGNE, "."), 
			inner.counts=NA, inner.iter=NA)
	res
}	


#pure fixed point iteration with a merit function
pureFP <- function(xinit, fn, merit, control, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	k <- 1
	xk_1 <- xinit
	xk <- fn(xk_1)
	xkp1 <- fn(xk)
	merit_xk <- merit(xk, y=xkp1)

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" m(x_k)", merit_xk,  "\n")
	
	
	while( abs( merit_xk ) > tol && k < maxit) 
	{
		xk_1 <- xk
		k <- k+1
	
		xk <- xkp1
		xkp1 <- fn(xk)
		merit_xk <- merit(xk, y=xkp1)
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" m(x_k)", merit_xk,  "\n")
	}
	
	list(par = xk, value=merit_xk , counts=c(fn=k+1, merit=k+1), iter = k, 
		 code=(k >= maxit)*4 + (abs(merit_xk) <= tol)*1 + (max(abs(xk - xk_1)/abs(xk)) <= tol)*2)
}


#Uryasev and Rubinstein (non optimized stepsize)
relaxationAlgoUR <- function(xinit, stepfunc, fn, merit, control, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
		
	
	k <- 1
	xk_1 <- xinit
	alphak <- stepfunc(k)
	xk <- (1-alphak) * xk_1 + alphak * fn(xk_1)
	if(!is.null(merit))
		merit_xk <- merit(xk)
	else
		merit_xk <- max(abs(xk - xk_1))
		
	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" m(x_k)", merit_xk,  "\n")
	if(echo >= 3)
		cat("step size", alphak, "\n\n")
	
	while( abs( merit_xk ) > tol && k < maxit)
	{
		xk_1 <- xk
		k <- k+1
		
		alphak <- stepfunc(k)
		xk <- (1-alphak) * xk_1 + alphak * fn(xk_1)
		if(!is.null(merit))
			merit_xk <- merit(xk)
		else
			merit_xk <- max(abs(xk - xk_1))
		
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" m(x_k)", merit_xk,  "\n")
		if(echo >= 3)
			cat("step size", alphak, "\n\n")
	}
	
	list(par = xk, value=merit_xk , counts=c(fn=k+1, merit=k+1), iter = k, 
		 code=(k >= maxit)*4 + (abs(merit_xk) <= tol)*1 + (max(abs(xk - xk_1)/abs(xk)) <= tol)*2)
}


#von Heusinger 
relaxationAlgoVH <- function(xinit, fn, merit, control, ...)
{
	sigma <- control$sigma
	beta <- control$beta
	tol <- control$tol
	echo <- control$echo
	maxit <- control$maxit
	
	k <- 0
	xk_1 <- xinit
	xk <- fn(xk_1)
	merit_xk <- merit(xk)
	
	counts <- c(fn=1, merit=1)

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" m(x_k)", merit_xk,  "\n")

	while( abs( merit_xk ) > tol && k < maxit)
	{
		k <- k+1
		xk_1 <- xk
		
		dk <- fn(xk) - xk
		normdk <- sqrt(sum(dk^2))
		merit_xk <- merit(xk)
		
		#backtracking line search
		tk <- 1
		l <- 0
		merit_xktkdk <- merit(xk + tk*dk)
		# use remark 4.32 
		while( merit_xktkdk > merit_xk - sigma * tk^2 * normdk^2 )
		{
			tk <- tk * beta
			merit_xktkdk <- merit(xk + tk*dk)
			if(echo >= 3)
			{
				cat(l, "\t", merit_xktkdk, "\t <= ")
				cat(merit_xk - sigma * tk^2 * normdk^2, "?\t")
				cat("tk", tk, "\n\n")
			}
			l <- l+1
		}
						
		xk <- xk_1 + tk * dk
		merit_xk <- merit_xktkdk
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" m(x_k)", merit_xk,  "\n")

		counts <- counts+1
		counts[2] <- counts[2]+l
	}	
	
	list(par = xk, value=merit_xk , counts=counts, iter = k, 
		 code=(k >= maxit)*4 + (abs(merit_xk) <= tol)*1 + (max(abs(xk - xk_1)/abs(xk)) <= tol)*2)
}



#extrapolation method for fixed point iteration
extrapolFP <- function(xinit, fn, merit, control, method, ...)
{
	echo <- control$echo
	maxit <- control$maxit
	tol <- control$tol
	
	k <- 1
	xk_1 <- xinit
	xk <- fn(xk_1)
	merit_xk <- merit(xk)
	
	counts <- c(fn=1, merit=1)

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" m(x_k)", merit_xk,  "\n")
	
	while( abs( merit_xk ) > tol && k < maxit) 
	{
		xk_1 <- xk
		k <- k+1
		
		xk <- fn(xk_1)
		xkp1 <- fn(xk)
		#RRE/MPE cycle of order 1, equivalent to Aitken Delta process when xk is univariate
		Delta1_xk <- xk - xk_1
		Delta2_xk <- xkp1 - 2*xk + xk_1

		if(method == "RRE" || method == "SqRRE")
			mystep <- crossprod(Delta2_xk, Delta1_xk) / crossprod(Delta2_xk, Delta2_xk)
		if(method == "MPE" || method == "SqMPE")
			mystep <- crossprod(Delta1_xk, Delta1_xk) / crossprod(Delta1_xk, Delta2_xk)

		if(is.nan(mystep) || is.infinite(mystep))
		{
			warning("Delta square xk is too small: Use the last iterate.")
			xk <- xkp1
		}else if(substr(method, 1, 2) == "Sq")
		{
			xk <- xk_1 - 2*Delta1_xk * mystep + Delta2_xk * mystep^2
		}else #simple version
		{
			xk <- xk_1 - Delta1_xk * mystep
		}
		
		merit_xk <- merit(xk)
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" m(x_k)", merit_xk,  "\n")
		if(echo >= 3)		
			cat(" Delta(xk)", Delta1_xk, "\n", "Delta^2(xk)", Delta2_xk, "\n\n")
		
		counts <- counts + c(2,1)
	}
	
	list(par = xk, value=merit_xk , counts=counts, iter = k, 
		 code=(k >= maxit)*4 + (abs(merit_xk) <= tol)*1 + (max(abs(xk - xk_1)/abs(xk)) <= tol)*2)
}


