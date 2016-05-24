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
### utility functions for evaluation with list-arguments in GNE
###
###         R functions
### 



evalwitharglist <- function(f, x, addarg)
{
	if(length(addarg) == 0)
	{	
		f(x)
	}else
	{	
		do.call(f, c(list(x), addarg) )
	}	
}


catlist <- function(x, ..., issummary=TRUE)
{
	if(issummary)
		cat(capture.output(summary(x)), sep="\n", ...)
	else
		cat(capture.output(print(x)), sep="\n", ...)
}



list2array <- function(x, listname=1:length(x))
array(unlist(x), 
	dim = c(dim(x[[1]]), length(x)),
	dimnames=c(dimnames(x[[1]]), list(listname))
) 

testfunc <- function(f, x=NULL, arg=NULL, echo=TRUE, errmess="", tobelisted=TRUE)
{
	nbtotarg <- length(formals(f))
	if(is.logical(echo))
		echo <- 1
	if(echo >= 2)
		print(formals(f))
	
	if(!is.null(arg) && tobelisted)
		arg <- list(arg)
	
	if(!is.null(x))
	{
		if(!is.null(arg))
		{
			if(nbtotarg < 2)
				stop("wrong number of arguments.")
			else if(nbtotarg == 2)
				arglist <- c(list(x), arg)
			else if(nbtotarg >= 3)
				arglist <- c(list(x), as.list(rep(1, nbtotarg-2)), arg)
		}else
		{
			if(nbtotarg < 1)
				stop("wrong number of arguments.")
			else if(nbtotarg == 1)
				arglist <- list(x)
			else if(nbtotarg >= 2)
				arglist <- c(list(x), as.list(rep(1, nbtotarg-1)))
		}
	}else
	{
		if(!is.null(arg))
		{
			if(nbtotarg < 1)
				stop("wrong number of arguments.")
			else if(nbtotarg == 1)
				arglist <- arg
			else if(nbtotarg >= 2)
				arglist <- c(as.list(rep(1, nbtotarg-1)), arg)
		}else
		{
			if(nbtotarg == 0)
				arglist <- list()
			else if(nbtotarg >= 1)
				arglist <- as.list(rep(1, nbtotarg))
		}		
	}
	if(echo >= 2)
	{
		print(is.function(f))
		print(f)
		print(arglist)
	}
	if(length(arglist) != nbtotarg)
		stop("wrong number of arguments when testing function.")
	
	test.try <- try(do.call(f, arglist), silent=echo >= 2)
					
	if(class(test.try) == "try-error")
	{
		if(echo >= 1)
		{
			cat("Error when calling function, below the try output.\n")
			print(test.try)
			cat("Arguments are:\n")
			print(arglist)
			stop(errmess)
		}
	}else if (echo >= 2)
		print(test.try)
	
	invisible()
}


