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
### utility functions for constrained equations in GNE
###
###         R functions
### 


#functions of the Constrained Equation Reformulation of the GNEP
#z = (x, lambda, mu, w, y) 
#or (x, lambda tilde, w tilde) with lambda tilde = (lambda, mu) and w tilde = (w, y)



testargfunCER <- function(z, dimx, dimlam, 
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	echo=FALSE)
{
	
	#sanity check	
	nplayer <- length(dimx)	
	if(!is.function(grobj))
		stop("argument grobj must be a function.")
	
	if(missing(constr))
		constr <- NULL
	else if(!is.function(constr) && !is.null(constr))
		stop("argument constr must be a function.")
	if(!is.null(constr) && missing(grconstr))
		stop("argument grconstr incompatible with constr.")
	else if(!is.null(constr) && !is.function(grconstr))
		stop("argument grconstr must be a function.")
	if(is.null(constr))
		dimlam <- rep(0, nplayer)
	else if(!is.null(constr) && (missing(dimlam) || length(dimlam) != nplayer))
		stop(paste("argument dimlam must be a vector of length", nplayer, "."))
	
	if(missing(joint))
		joint <- NULL
	else if(!is.function(joint) && !is.null(joint))
		stop("argument joint must be a function.")
	if(!is.null(joint) && missing(grjoint))
		stop("argument grjoint incompatible with joint.")
	else if(!is.null(joint) && !is.function(grjoint))
		stop("argument grjoint must be a function.")
	if(is.null(joint))
		dimmu <- 0
	else if(!is.null(joint) && (missing(dimmu) || length(dimmu) != 1))
		stop("argument dimmu must be a vector of length 1.")
	
	if(echo)
	{
		print(dimx)
		print(dimlam)
		print(length(z))
	}
	
	if(length(z) != sum(dimx) + 2*sum(dimlam) + 2*sum(dimmu))
		stop("CER: incompatible dimension for dimlam, dimx, dimmu.")		
		
	#objective gradient
	if(!missing(arggrobj) && !is.null(arggrobj))
		grobjfinal <- grobj
	else
	{
		grobjfinal <- function(z, i, j, ...) grobj(z, i, j)
		arggrobj <- list()				   
	}	
	str <- paste("the call to grobj(z, 1, 1, arggrobj) does not work.", "arguments are", 
				 paste(names(formals(grobj)), collapse=","), ".")
	testfunc(grobjfinal, z, arg=arggrobj, echo=echo, errmess=str)
	
	
	#constraint	
	if(!is.null(constr))
	{
		if(!missing(argconstr) && !is.null(argconstr))
			constrfinal <- constr
		else
		{
			constrfinal <- function(z, i, ...) constr(z, i)
			argconstr <- list()	
		}
		str <- paste("the call to constr(z, 1, argconstr) does not work.", "arguments are", 
					 paste(names(formals(constr)), collapse=","), ".")
		testfunc(constrfinal, z, arg=argconstr, echo=echo, errmess=str)
		
		#constraint gradient
		if(!missing(arggrconstr) && !is.null(arggrconstr))
			grconstrfinal <- grconstr
		else
		{
			grconstrfinal <- function(z, i, j, ...) grconstr(z, i, j)
			arggrconstr <- list()
		}
		str <- paste("the call to grconstr(z, 1, 1, arggrconstr) does not work.", "arguments are", 
					 paste(names(formals(grconstr)), collapse=","), ".")
		testfunc(grconstrfinal, z, arg=arggrconstr, echo=echo, errmess=str)
	}else
	{
		constrfinal <- grconstrfinal <- argconstr <- arggrconstr <- NULL
	}
		
	
	#joint constraint function
	if(!is.null(joint))
	{
		if(!missing(argjoint) && !is.null(argjoint))
			jointfinal <- joint
		else
		{
			jointfinal <- function(z, ...) joint(z)
			argjoint <- list()	
		}
		str <- paste("the call to joint(z, argconstr) does not work.", "arguments are", 
					 paste(names(formals(joint)), collapse=","), ".")
		testfunc(jointfinal, z, arg=argjoint, echo=echo, errmess=str)
		
		#gradient
		if(!missing(arggrjoint) && !is.null(arggrjoint))
			grjointfinal <- grjoint
		else
		{
			grjointfinal <- function(z, j, ...) grjoint(z, j)
			arggrjoint <- list()
		}
		str <- paste("the call to grjoint(z, 1, arggrjoint) does not work.", "arguments are", 
					 paste(names(formals(grjoint)), collapse=","), ".")
		testfunc(grjointfinal, z, arg=arggrjoint, echo=echo, errmess=str)
	}else
	{
		jointfinal <- grjointfinal <- argjoint <- arggrjoint <- NULL
	}
	
	list(grobj=grobjfinal, arggrobj=arggrobj, constr=constrfinal, argconstr=argconstr,
		 grconstr=grconstrfinal, arggrconstr=arggrconstr, joint=jointfinal, 
		 argjoint=argjoint, grjoint=grjointfinal, arggrjoint=arggrjoint,
		 dimx=dimx, dimlam=dimlam, dimmu=dimmu, nplayer=nplayer)
}

testargjacCER <- function(z, dimx, dimlam,
	heobj, argheobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	hejoint, arghejoint,
	echo=FALSE)
{
	
	#sanity check	
	nplayer <- length(dimx)	
	if(!is.function(heobj))
		stop("argument heobj must be a function.")
	
	if(missing(constr))
		constr <- NULL
	else if(!is.function(constr) && !is.null(constr))
		stop("argument constr must be a function.")
	if(!is.null(constr) && missing(grconstr))
		stop("argument grconstr incompatible with constr.")
	else if(!is.null(constr) && !is.function(grconstr))
		stop("argument grconstr must be a function.")
	else if(!is.null(constr) && !is.function(heconstr))
		stop("argument heconstr must be a function.")
	if(is.null(constr))
		dimlam <- rep(0, nplayer)
	else if(!is.null(constr) && (missing(dimlam) || length(dimlam) != nplayer))
		stop(paste("argument dimlam must be a vector of length", nplayer, "."))
	
	if(missing(joint))
		joint <- NULL
	else if(!is.function(joint) && !is.null(joint))
		stop("argument joint must be a function.")
	if(!is.null(joint) && missing(grjoint))
		stop("argument grjoint incompatible with joint.")
	else if(!is.null(joint) && !is.function(grjoint))
		stop("argument grjoint must be a function.")
	else if(!is.null(joint) && !is.function(hejoint))
		stop("argument hejoint must be a function.")
	if(is.null(joint))
		dimmu <- 0
	else if(!is.null(joint) && (missing(dimmu) || length(dimmu) != 1))
		stop("argument dimmu must be a vector of length 1.")
	
	if(length(z) != sum(dimx) + 2*sum(dimlam) + 2*sum(dimmu))
		stop("incompatible dimension for dimlam, dimx, dimmu.")		
	
	
	#objective Hessian
	if(!missing(argheobj) && !is.null(argheobj))
		heobjfinal <- heobj
	else
	{
		heobjfinal <- function(z, i, j, k, ...) heobj(z, i, j, k)
		argheobj <- list()				   
	}	
	str <- paste("the call to heobj(z, 1, 1, 1, argheobj) does not work.", "arguments are", 
				 paste(names(formals(heobj)), collapse=","), ".")
	testfunc(heobjfinal, z, arg=argheobj, echo=echo, errmess=str)
	
	#constraint	
	if(!is.null(constr))
	{
		if(!missing(argconstr) && !is.null(argconstr))
		constrfinal <- constr
		else
		{
			constrfinal <- function(z, i, ...) constr(z, i)
			argconstr <- list()	
		}
		str <- paste("the call to constr(z, 1, argconstr) does not work.", "arguments are", 
					 paste(names(formals(constr)), collapse=","), ".")
		testfunc(constrfinal, z, arg=argconstr, echo=echo, errmess=str)
		
		#constraint gradient
		if(!missing(arggrconstr) && !is.null(arggrconstr))
		grconstrfinal <- grconstr
		else
		{
			grconstrfinal <- function(z, i, j, ...) grconstr(z, i, j)
			arggrconstr <- list()
		}
		str <- paste("the call to grconstr(z, 1, 1, arggrconstr) does not work.", "arguments are", 
					 paste(names(formals(grconstr)), collapse=","), ".")
		testfunc(grconstrfinal, z, arg=arggrconstr, echo=echo, errmess=str)

		#constraint hessian
		if(!missing(argheconstr) && !is.null(argheconstr))
			heconstrfinal <- heconstr
		else
		{
			heconstrfinal <- function(z, i, j, k, ...) heconstr(z, i, j, k)
			argheconstr <- list()
		}
		str <- paste("the call to heconstr(z, 1, 1, 1, argheconstr) does not work.", "arguments are", 
					 paste(names(formals(heconstr)), collapse=","), ".")
		testfunc(heconstrfinal, z, arg=argheconstr, echo=echo, errmess=str)
	}else
	{
		constrfinal <- grconstrfinal <- heconstrfinal <- argconstr <- arggrconstr <- argheconstr <- NULL
	}	
	
	#joint constraint function
	if(!is.null(joint))
	{
		if(!missing(argjoint) && !is.null(argjoint))
			jointfinal <- joint
		else
		{
			jointfinal <- function(z, ...) joint(z)
			argjoint <- list()	
		}
		str <- paste("the call to joint(z, argconstr) does not work.", "arguments are", 
					 paste(names(formals(joint)), collapse=","), ".")
		testfunc(jointfinal, z, arg=argjoint, echo=echo, errmess=str)
		
		#gradient
		if(!missing(arggrjoint) && !is.null(arggrjoint))
			grjointfinal <- grjoint
		else
		{
			grjointfinal <- function(z, j, ...) grjoint(z, j)
			arggrjoint <- list()
		}
		str <- paste("the call to grjoint(z, 1, arggrjoint) does not work.", "arguments are", 
					 paste(names(formals(grjoint)), collapse=","), ".")
		testfunc(grjointfinal, z, arg=arggrjoint, echo=echo, errmess=str)
		
		#hessian
		if(!missing(arghejoint) && !is.null(arghejoint))
			hejointfinal <- hejoint
		else
		{
			hejointfinal <- function(z, j, k, ...) hejoint(z, j, k)
			arghejoint <- list()
		}
		str <- paste("the call to hejoint(z, 1, 1, arghejoint) does not work.", "arguments are", 
					 paste(names(formals(hejoint)), collapse=","), ".")
		testfunc(hejointfinal, z, arg=arghejoint, echo=echo, errmess=str)

	}else
	{
		jointfinal <- grjointfinal <- hejointfinal <- argjoint <- arggrjoint <- arghejoint <- NULL
	}
	
	list(heobj=heobjfinal, argheobj=argheobj, constr=constrfinal, argconstr=argconstr,
		 grconstr=grconstrfinal, arggrconstr=arggrconstr, heconstr=heconstrfinal, 
		 argheconstr=argheconstr, joint=jointfinal, argjoint=argjoint, grjoint=grjointfinal, 
		 arggrjoint=arggrjoint, hejoint=hejointfinal, arghejoint=arghejoint,
		 dimx=dimx, dimlam=dimlam, dimmu=dimmu, nplayer=nplayer)	
}

