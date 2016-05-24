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
### utility functions for rejection algorithm in GNE
###
###         R functions
### 



#rejection algorithm with different distributions

rejection <- function(constr, nvars, LB=0, UB=1, ..., echo=FALSE, 
	method=c("unif","norm", "normcap"), control=list())
{
	#default control parameters
	con <- list(mean=(LB + UB)/2, sd=-(UB-LB) /4 /qnorm(0.025))
	namc <- names(con)
	con[namc <- names(control)] <- control
	
	
	#adjust the lentgh of lower/upper bounds
	LB <- rep(LB, length=nvars)	
	method <- match.arg(method, c("unif","norm","normcap"))
	
	if(method == "unif")
	{
		#draw
		x <- runif(nvars, LB, UB)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) )
		{
			if(echo) print(x)
			x <- runif(nvars, LB, UB)
		}
		return(x)
	}
	
	if(method == "norm")
	{
		#draw
		x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) )
		{
			if(echo) print(x)
			x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		}
		return(x)
	}
	
	if(method == "normcap")
	{
		#draw
		x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		
		#until all constraints are satisfied
		while( any( constr(x, ...) >= 0 ) || sum( (UB-x)^2 ) <= 1)
		{
			if(echo) print(x)
			x <- rnorm(nvars, mean=con$mean, sd=con$sd)
		}
		return(x)
	}
	
	stop("wrong method for rejection algo.")
	
}

	