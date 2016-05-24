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
### utility functions for minimization in GNE
###
###         R functions
### 



#steepest descent methods
minfinding <- function(xinit, func, gradfunc, method, control, ...)
{
	tol <- control$tol
	maxit <- control$maxit
	echo <- control$echo
	stepinit <- control$stepinit
	k <- 0		
	
	if(missing(func) || missing(xinit) || missing(gradfunc))
		stop("Missing parameters.")
	
	testmyfunc <- func(xinit, ...)
	testmygrad <- gradfunc(xinit, ...)
	
	if(class(testmyfunc) != class(testmygrad))
		stop("Wrong class for argument results.")
	
	
	if(class(testmyfunc) %in% c("integer", "numeric"))
	{	
		inner.counts.fn <- inner.counts.gr <- NULL
		inner.iter.fn <- inner.iter.gr <- NULL
		noitercount <- TRUE
	}else if(class(testmyfunc) == "list")
	{
		noitercount <- FALSE
		
		inner.counts.fn <- c(0, 0)
		inner.iter.fn <- 0
		wrapfunc <- function(x, ...) 
		{
			fx <- func(x, ...)
			inner.counts.fn <<- inner.counts.fn + fx$counts
			inner.iter.fn <<- inner.iter.fn + fx$iter
			fx$value
		}
		
		inner.counts.gr <- c(0, 0)
		inner.iter.gr <- 0
		wrapgrad <- function(x, ...) 
		{
			gx <- gradfunc(x, ...)
			inner.counts.gr <<- inner.counts.gr + gx$counts
			inner.iter.gr <<- inner.iter.gr + gx$iter
			gx$value
		}
	}else
		stop("Unsupported class for argument functions.")
	
	if(method == "BFGS")
	{
		if(noitercount)
			resoptim <- optim(xinit, func, gradfunc, method="BFGS", 
							  control=list(maxit=maxit, trace=echo, REPORT=1, abstol=tol), ...)
		else
			resoptim <- optim(xinit, wrapfunc, wrapgrad, method="BFGS", 
							  control=list(maxit=maxit, trace=echo, REPORT=1, abstol=tol), ...)
		
		if(echo >=2)
		{
			print(resoptim)
			cat("\n")
		}	
		k <- as.numeric(resoptim$counts[2])-1
		xk <- as.numeric(resoptim$par)
		f_xk <- as.numeric(resoptim$value)
		
		counts <- resoptim$counts
		names(counts) <- c("Vhat", "gradVhat")

	}	
	
	if(method == "BB")
	{		
		if(noitercount)
			resBB <- BBmethod(xinit, func, gradfunc, control=control, ...)
		else	
			resBB <- BBmethod(xinit, wrapfunc, wrapgrad, control=control, ...)
		
		k <- as.numeric(resBB$iter)-1
		xk <- as.numeric(resBB$par)
		f_xk <- as.numeric(resBB$value)

		counts <- resBB$counts
	}
	
	
	list(par=xk, outer.counts=counts, outer.iter=k, code=(k >= maxit)*1 + (abs(f_xk) > tol)*10, 
		 inner.iter.fn=inner.iter.fn, inner.iter.gr=inner.iter.gr, inner.counts.fn=inner.counts.fn, 
		 inner.counts.gr=inner.counts.gr)
}


BBmethod <- function(xinit, func, gradfunc, control, ...)
{
	tol <- control$tol
	maxit <- control$maxit
	echo <- control$echo
	stepinit <- control$stepinit
	
	k <- 0		
	xk_1 <- xinit
	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk_1, "\n")
	
	k <- k+1
	xk <- xk_1 - stepinit * gradfunc(xk_1, ...)
	f_xk <- func(xk, ...)

	if(echo >= 1)
		cat("**** k", k, "\n x_k", xk, "\n")
	if(echo >= 2)
		cat(" ||V(x_k)||", f_xk,  "\n")
	
	while(abs(f_xk) > tol && k < maxit)		
	{
		k <- k+1
		xkp1 <- BBnext(gradfunc, xk, xk_1, ...)
		xk_1 <- xk
		xk <- xkp1
		f_xk <- func(xk, ...)
		
		if(echo >= 1)
			cat("**** k", k, "\n x_k", xk, "\n")
		if(echo >= 2)
			cat(" ||V(x_k)||", f_xk,  "\n")
		
	}
	counts <- c(1,1) + c(k, 2*k)
	names(counts) <- c("Vhat", "gradVhat")
	
	list(par=xk, value=f_xk, counts=counts, iter=k, code=(k >= maxit)*1 + (abs(f_xk) > tol)*10)
}


#Barzilai Borwein
BBnext <- function(grad, xk, xk_1, ...)
{
	diff <- xk - xk_1	
	gk <- grad(xk, ...)
	step <- gk - grad(xk_1, ...)

	#equation 2.3 of Fletcher (2001)
	slope <- as.numeric( (t(step) %*% diff) / (t(step) %*% step) )
	
#	#equation 2.2 of Fletcher (2001)
#	slope <- as.numeric( (t(diff) %*% step) / (t(diff) %*% diff) )
	
	xk - slope * gk
}


#Newton step
NewtonNext <- function(xk, f_xk, jacf_xk, silent=TRUE)
{
	A <- jacf_xk
	b <- -f_xk	
	mycatch <- try( direction <- qr.solve(A, b) , silent=silent)
	
	if(class(mycatch) == "try-error")
	{	
		if(!silent)
			cat(mycatch)
		stop("The Newton step cannot be done.")
	}else
		xk + direction
}


#Levenberg-Marquardt step
LevenMarqNext <- function(xk, f_xk, jacf_xk, silent=TRUE, delta=1)
{
	lambdak <- sqrt( sum( f_xk^2 ) )^delta
	A <- crossprod( jacf_xk, jacf_xk ) + lambdak * diag( length(xk) )
	b <- -crossprod( jacf_xk, f_xk )
	mycatch <- try( direction <- solve(A, b) , silent=silent)
	
	if(class(mycatch) == "try-error")
	{	
		if(!silent)
			cat(mycatch)
		stop("The Levenberg-Marquart step cannot be done.")
	}else
		xk + as.vector( direction )
}

#Quasi Newton step
QuasiNewtonNext <- function(xk, f_xk, Wk, silent=TRUE, inv=TRUE)
{
	if(inv)
		mycatch <- try( direction <- - Wk %*% f_xk , silent=silent)
	if(!inv)
		mycatch <- try( direction <- solve(Wk, -f_xk), silent=silent)
		
	if(class(mycatch) == "try-error")
	{	
		if(!silent)
			cat(mycatch)
		stop("The quasi-Newton step cannot be done.")
	}else
		xk + as.vector( direction )
}
