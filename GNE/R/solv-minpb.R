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
### minimization in GNE
###
###         R functions
### 


#algorithm for minimum finding of function fn, with gradient gr.
#possibly with a constraint function hin(x)<=0 with Jacobian hin.jac
minpb <- function(xinit, fn, argfn, gr=NULL, arggr, 
	hin=NULL, arghin, hin.jac=NULL, arghin.jac, method=c("BB","BFGS", "CG", "Brent"), 
	control=list(), silent=TRUE, ...)
{
	method <- match.arg(method)
	if(missing(argfn))
		argfn <- NULL
	if(missing(arggr))
		arggr <- NULL
	if(missing(arghin))
		arghin <- NULL
	if(missing(arghin.jac))
		arghin.jac <- NULL
	
	noitercount <- FALSE
	#inner iterations to compute minimization function fn
	inner.counts.fn <- c(0, 0) #call to gap and grad gap
	inner.iter.fn <- 0 #iter number
	
	finit <- evalwitharglist(fn, xinit, argfn)
	
	if(class(finit) == "list")
	Fn <- function(x) 
	{
		fx <- evalwitharglist(fn, x, argfn)
		inner.counts.fn <<- inner.counts.fn + fx$counts
		inner.iter.fn <<- inner.iter.fn + fx$iter
		fx$value
	}
	else
	Fn <- function(x) 
		evalwitharglist(fn, x, argfn)
	
	
	
	#default control parameters
	conBB <- list(ftol=1e-8, maxit=100, trace=0, gtol = 1e-05, maxfeval = 10000, M=5,
				  maximize = FALSE, triter = 1, eps = 1e-07, checkGrad=FALSE)
	if(!is.null(control$tol))
		conBB$ftol <- control$tol
	if(!is.null(control$maxit))
		conBB$maxit <- control$maxit
	if(!is.null(control$trace))
		conBB$trace <- control$trace
	if(!is.null(control$checkGrad))
		conBB$checkGrad <- control$checkGrad	
	
	conoptim <- list(abstol=1e-8, maxit=100, trace=0, REPORT=1)
	if(!is.null(control$tol))
		conoptim$abstol <- control$tol
	if(!is.null(control$maxit))
		conoptim$maxit <- control$maxit
	if(!is.null(control$trace))
		conoptim$trace <- control$trace
	if(!is.null(control$REPORT))
		conoptim$REPORT <- control$REPORT
	conalabama <- list(eps=1e-8, itmax=100, trace=FALSE)
	if(!is.null(control$tol))
		conalabama$eps <- control$tol
	if(!is.null(control$maxit))
		conalabama$itmax <- control$maxit
	if(!is.null(control$trace))
		conalabama$trace <- control$trace
	
	if(missing(fn) || missing(xinit))
		stop("missing parameters.")

	#create the objective function
	test.try <- try( Fn(xinit), silent=silent )
	if(class(test.try) == "try-error")
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="Can't evalate fn(init).", fvec=NA) )
	if(is.na(test.try) || unlist(sapply(test.try, is.nan)) || unlist(sapply(test.try, is.infinite)))
		return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				 message="fn(init) has infinite, NA or NaN values.", fvec=NA) )
	if(!silent)
		cat("test of function done.\n")
		
	
	#create the gradient function (if provided)
	if(!is.null(gr))
	{
		grinit <- evalwitharglist(gr, xinit, arggr)
		
		#inner iterations to compute gradient function gr
		inner.counts.gr <- c(0, 0) #call to gap and grad gap
		inner.iter.gr <- 0 #iter number
		if(class(grinit) == "list")
			Gr <- function(x)
			{
				g <- evalwitharglist(gr, x, arggr)
				inner.counts.gr <<- inner.counts.gr + g$counts
				inner.iter.gr <<- inner.iter.gr + g$iter
				g$value
			}
		else
			Gr <- function(x) evalwitharglist(gr, x, arggr)
	}else
	{
		Gr <- NULL
		inner.counts.gr <- inner.iter.gr <- 0
	}
	
	#create the inequality constraint g(x) >= 0 (to pass constrOptim.nl) while in argument joint(x) <= 0
	if(!is.null(hin))
	{
		Hin <- function(x) -evalwitharglist(hin, x, arghin)
		test.try <- try( Hin(xinit), silent=silent )
		if(class(test.try) == "try-error")
		{	
			Hin <- function(x) -evalwitharglist(hin, x, list(arghin))
			test.try <- try( Hin(xinit), silent=silent )			
		}
		if(class(test.try) == "try-error")
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="Can't evalate hin(init).", fvec=NA) )
		if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="hin(init) has infinite, NA or NaN values.", fvec=NA) )
	}else
		Hin <- NULL
	if(!silent)
		cat("test of inequality constraint done.\n")
	if(method == "BB" && !is.null(Hin))
		stop("BB does not allow constraints.")
	
	
	#create the inequality constraint Jacobian
	if(!is.null(hin.jac))
	{
		Hin.jac <- function(x) -evalwitharglist(hin.jac, x, arghin.jac)
		test.try <- try( Hin.jac(xinit), silent=silent )
		if(class(test.try) == "try-error")
		{	
			Hin.jac <- function(x) -evalwitharglist(hin.jac, x, list(arghin.jac))
			test.try <- try( Hin.jac(xinit), silent=silent )			
		}
		if(class(test.try) == "try-error")
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="Can't evalate jac.hin(init).", fvec=NA) )
		if(any(is.na(test.try) || is.nan(test.try) || is.infinite(test.try)) )
			return( list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					 message="jac.hin(init) has infinite, NA or NaN values.", fvec=NA) )		
	}else
		Hin.jac <- NULL

	

	if(!silent)
		cat("beginning of optimization:\t")
	
	if(is.null(Hin))
	{
		if(method != "BB")
		{
			if(!silent)
				cat("optimization carried out by optim\n")
			test.try <- try( optim(par=xinit, fn=Fn, gr=Gr, control=conoptim, 
							   method=method, ...), silent=silent)
		}else
		{
			if(!silent)
				cat("optimization carried out by BBoptim\n")
			test.try <- try( spg(par=xinit, fn=Fn, gr=Gr, control=conBB, method=1, ...), silent=silent)
		}
	}else 
	{
		
		if(!silent)
			cat("optimization carried out by constrOptim.nl\n")
		test.try <- try( constrOptim.nl(par=xinit, fn=Fn, gr=Gr, hin=Hin, 
						hin.jac=Hin.jac, control.outer=conalabama, ...), silent=silent)
		if(class(test.try) == "try-error")
			res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
					message=paste("Error in the non smooth problem:", test.try, "."))
		
	}
	if(!silent)
		cat("end of optimization\n")
	
	if(class(test.try) == "try-error")
	{
		res <- list(par= NA, value=NA, counts=NA, iter=NA, code=100, 
				message=paste("Error in the non smooth problem:", test.try, "."))
	}else
	{
		if(method != "BB" && is.null(Hin))
		{
			code <- ifelse(test.try$convergence == 0, 1, ifelse(test.try$convergence == 1, 4, 2))
			res <- list(par = test.try$par, value = test.try$value, 
				outer.counts = test.try$counts, outer.iter = as.numeric(test.try$counts["function"]), 
				inner.counts = inner.counts.fn+inner.counts.gr, inner.iter = inner.iter.fn+inner.iter.gr, 
				code = code, message = test.try$message)		
		}
		if(method == "BB" && is.null(Hin))
		{
			code <- ifelse(test.try$convergence == 0, 1, ifelse(test.try$convergence == 1, 4, 2))
			res <- list(par = test.try$par, value = test.try$value, 
				outer.counts = c(test.try$feval, test.try$iter+1), outer.iter = test.try$iter, 
				inner.counts = inner.counts.fn+inner.counts.gr, inner.iter = inner.iter.fn+inner.iter.gr, 
				code = code, message = test.try$message)					
		}
	
		if(!is.null(Hin))
		{
			code <- ifelse(test.try$convergence == 0, 1, ifelse(test.try$convergence == 7, 4, 2))
			res <- list(par = test.try$par, value = test.try$value, 
				outer.counts = test.try$counts, outer.iter = test.try$outer.iterations, 
				inner.counts = inner.counts.fn+inner.counts.gr, inner.iter = inner.iter.fn+inner.iter.gr, 
				code = code, message = test.try$message)		
		}
	}
	
	res
}



