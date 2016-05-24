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
### utility functions for variational inequality Reformulation in GNE
###
###         R functions
### 


#functions of the Variational Inequality Reformulation of the GNEP


gapVIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	#sanity check		
	arg <- testarggapVIR(x, dimx, grobj, arggrobj)
	
	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")

	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	
	#ith increment
	grobjregi <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		norm_xy <- sum( (x[idx_i] - y[idx_i])^2 )
		arg$grobj(x, i, i, arg$arggrobj) * (x[idx_i] - y[idx_i]) - par$alpha/2*norm_xy
	}
	sum( sapply(1:nplayer, grobjregi) )
}


gradxgapVIR <- function(x, y, dimx, grobj, arggrobj, 
	heobj, argheobj, param=list(), echo=FALSE)
{
	arg <- testarggradxgapVIR(x, dimx, grobj, arggrobj, 
		heobj, argheobj, echo)

	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")
	
	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	HeGr_i <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		GrObj_i <- sapply(idx_i, function(j) arg$grobj(x, i, j, arg$arggrobj))
		HessObj_i <- sapply(idx_i, function(k) sapply(idx_i, function(j) arg$heobj(x, i, j, k, arg$argheobj)))
		HessObj_i %*% (x[idx_i] - y[idx_i]) + GrObj_i - par$alpha * (x[idx_i] - y[idx_i])
	}
	sapply(1:nplayer, HeGr_i)
}

gradygapVIR <- function(x, y, dimx, grobj, arggrobj, param=list(), echo=FALSE)
{
	arg <- testarggradygapVIR(x, dimx, grobj, arggrobj, echo)
	
	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer
	
	if(length(x) != n)
		stop("wrong argument x")
	if(length(y) != n)
		stop("wrong argument y")
	
	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	Gry_i <- function(i)
	{
		idx_i <- index4x[1,i]:index4x[2,i]
		GrObj_i <- sapply(idx_i, function(j) arg$grobj(x, i, j, arg$arggrobj))
		-GrObj_i + par$alpha * (x[idx_i] - y[idx_i])
	}
	
	sapply(1:nplayer, Gry_i)
}

fpVIR <- function(x, dimx, obj, argobj, joint, argjoint,  
	grobj, arggrobj, jacjoint, argjacjoint, param=list(), 
	echo=FALSE, control=list(), yinit=NULL, optim.method="default")
{
	arg <- testargfpVIR(x, dimx, obj, argobj, joint, argjoint,  
			grobj, arggrobj, jacjoint, argjacjoint, echo)

	#default parameters
	par <- list(alpha = 1e-1)
	namp <- names(par)
	par[namp <- names(param)] <- param
	
	if(optim.method == "default")
		optim.method <- "BFGS"
	#default control parameters
	con1.nl <- list(trace = echo, eps=1e-6, method=optim.method, itmax=100, NMinit=FALSE)
	namc1 <- names(con1.nl)
	con1.nl[namc1 <- names(control)] <- control
	con2.un <- list(trace = echo, fnscale = 1, parscale = rep.int(1, length(x)),
		ndeps = rep.int(1e-3, length(x)), maxit = 100L, abstol = -Inf, 
		reltol = 1e-6, alpha = 1.0, beta = 0.5, gamma = 2.0, REPORT = 1,
		type = 1, lmm = 5, factr = 1e7, pgtol = 0, tmax = 10, temp = 10.0)
	namc2 <- names(con2.un)
	con2.un[namc2 <- names(control)] <- control
	
	dimx <- arg$dimx
	n <- sum(arg$dimx)
	nplayer <- arg$nplayer

	if(is.null(yinit) && !is.null(arg$joint))
	{
		yinit <- rejection(function(x) arg$joint(x, arg$argjoint), n, method="norm")
		if(echo)
			cat("initial point in fpVIR", yinit, "-", arg$joint(yinit, arg$argjoint), "\n")
	}else if(is.null(yinit) && is.null(arg$joint))
	{
		yinit <- rnorm(n)
	}else
	{
		yinit <- yinit[1:n]
	}
	
	if(is.null(arg$joint) && !is.null(arg$grobj))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param)
			gr <- function(y, z, param) -gradygapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param) 
			
			res <- optim(yinit, fn=fn, gr=gr, method=optim.method, 
						 z=x, param=par, control=con2.un)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.\n")
			c(res, optim.function="optim", optim.method=optim.method)
		}
	}else if(is.null(arg$joint) && is.null(arg$grobj))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param)
			
			res <- optim(yinit, fn=fn, method=optim.method, 
						 z=x, param=par, control=con2.un)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.\n")
			c(res, optim.function="optim", optim.method=optim.method)
		}
	}else if(!is.null(arg$joint) && !is.null(arg$jacjoint))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param)
			gr <- function(y, z, param) -gradygapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param) 
			hin <- function(y, z, param) -arg$joint(y, arg$argjoint)
			hin.jac <- function(y, z, param) -arg$jacjoint(y, arg$argjacjoint)
			
			res <- constrOptim.nl(yinit, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac, 
								  z=x, param=par, control.outer=con1.nl)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.\n")
#			if(echo)
#				print(res)
			c(res, optim.function="constrOptim.nl", optim.method=optim.method)
		}
	}else if(!is.null(arg$joint) && is.null(arg$jacjoint))
	{
		yofx <- function(x)
		{
			fn <- function(y, z, param) -gapVIR(z, y, dimx, arg$grobj, arg$arggrobj, param)
			hin <- function(y, z, param) -arg$joint(y, arg$argjoint)
			
			res <- constrOptim.nl(yinit, fn=fn, hin=hin, 
								  z=x, param=par, control.outer=con1.nl)
			if(echo && res$convergence != 0)
				cat("Unsuccessful convergence.\n")
			c(res, optim.function="constrOptim.nl", optim.method=optim.method)
		}
	}else
		stop("missing argument in obj, grobj, joint, jacjoint.")
	
	yofx(x)
}

