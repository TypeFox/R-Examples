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
### utility functions for Nikaido-Isoda reformulation in GNE
###
###         R functions
### 


#functions of the Nikaido Isoda Reformulation of the GNEP
#z = (x) 



testarggapNIR <- function(z, dimx, 
	obj, argobj,
	echo=FALSE)
{
	
	#sanity check	
	nplayer <- length(dimx)	
	if(!is.function(obj))
		stop("argument obj must be a function.")
			
	if(echo)
	{
		print(dimx)
		print(length(z))
	}
	
	if(length(z) != sum(dimx))
		stop("CER: incompatible dimension for dimx.")		
		
	#objective function
	if(!missing(argobj) && !is.null(argobj))
		objfinal <- obj
	else
	{
		objfinal <- function(z, i, ...) obj(z, i)
		argobj <- list()				   
	}	
	str <- paste("the call to obj(z, 1, argobj) does not work.", "arguments are", 
				 paste(names(formals(obj)), collapse=","), ".")
	testfunc(objfinal, z, arg=argobj, echo=echo, errmess=str)
	
	list(obj=objfinal, argobj=argobj, dimx=dimx, nplayer=nplayer)
}


testarggradgapNIR <- function(z, dimx, 
	grobj, arggrobj, 
	echo=FALSE)
{
	
	#sanity check	
	nplayer <- length(dimx)	
	if(!is.function(grobj))
		stop("argument grobj must be a function.")
	

	if(echo)
	{
		print(dimx)
		print(length(z))
	}
	
	if(length(z) != sum(dimx))
	stop("CER: incompatible dimension for dimx.")		

	
	
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
	
	
	list(grobj=grobjfinal, arggrobj=arggrobj, dimx=dimx, nplayer=nplayer)
}


testargfpNIR <- function(z, dimx, 
	obj, argobj, 
	joint, argjoint,  
	grobj, arggrobj,
	jacjoint, argjacjoint,
	echo=FALSE)
{
	
	#sanity check	
	nplayer <- length(dimx)		
	if(!is.function(obj))
		stop("argument obj must be a function.")

	if(missing(joint))
		joint <- NULL
	else if(!is.function(joint) && !is.null(joint))
		stop("argument joint must be a function.")

	if(length(z) != sum(dimx))
		stop("incompatible dimension for dimx.")		
	
	if(missing(grobj))
		grobj <- NULL
	if(missing(jacjoint))
		jacjoint <- NULL
#	if(is.null(jacjoint) && !is.null(grobj))
#		stop("missing jacjoint argument.")
	if(!is.null(jacjoint) && is.null(grobj))
		stop("missing grobj argument.")
	
	#objective function
	if(!missing(argobj) && !is.null(argobj))
		objfinal <- obj
	else
	{
		objfinal <- function(z, i, ...) obj(z, i)
		argobj <- list()				   
	}	
	str <- paste("the call to obj(z, 1, argobj) does not work.", "arguments are", 
				 paste(names(formals(obj)), collapse=","), ".")
	testfunc(objfinal, z, arg=argobj, echo=echo, errmess=str)

	
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
	}else
		jointfinal <- argjoint <- NULL
	
	if(!is.null(grobj))
	{
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
		
	}else
		grobjfinal <- arggrobj <- NULL
	
	if(!is.null(jacjoint))
	{
		if(!missing(argjacjoint) && !is.null(argjacjoint))
			jacjointfinal <- jacjoint
		else
		{
			jacjointfinal <- function(z, ...) jacjoint(z)
			argjacjoint <- list()	
		}
		str <- paste("the call to jacjoint(z, argjacjoint) does not work.", "arguments are", 
					 paste(names(formals(jacjoint)), collapse=","), ".")
		testfunc(jacjointfinal, z, arg=argjacjoint, echo=echo, errmess=str)				
	}else
		jacjointfinal <- argjacjoint <- NULL
	
	list(obj=objfinal, argobj=argobj, 
		 grobj=grobjfinal, arggrobj=arggrobj, 
		 joint=jointfinal, argjoint=argjoint, 
		 jacjoint=jacjointfinal, argjacjoint=argjacjoint, 
		 dimx=dimx, nplayer=nplayer)
}

