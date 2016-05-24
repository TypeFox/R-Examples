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
### utility functions for line-search techniques in GNE
###
###         R functions
### 



# geometric backtracking line search - minstep based on Dennis and Schnabel
linesearch.geom <- function(xk, dk, slopek, con, merit, ...)
{
	stepk <- 1
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
	
	normfp <- merit(xk + stepk*dk, ...)
	normfk <- merit(xk, ...)
			
	while(stepk > minstep)
	{
		normfp <- merit(xk + stepk*dk, ...)
		
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
		
		#check Armijo condition
		if(normfp <= normfk + con$btol * stepk * slopek)
		{
			break
		}
		
		inner.iter <- inner.iter + 1	
		stepk <- con$sigma * stepk
	}
	if(con$trace >= 3)
		cat("xk", xk, "\ndk", dk, "stepk", stepk, "normfk", normfk, "normfp", normfp, "\n")
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)

}			

# geometric backtracking line search customized for interior point methods
linesearch.geom.cond <- function(xk, dk, slopek, con, merit, checkint, dimx, dimlam, ...)
{
	stepk <- 1
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		

	normfk <- merit(xk, dimx, dimlam, ...)	
	
	while(stepk > minstep)
	{
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\n")			
		if(checkint(xk + stepk*dk, dimx, dimlam))
		{	
			normfp <- merit(xk + stepk*dk, dimx, dimlam, ...)
			
			#traces in R console	
			if(con$trace >= 3)
				cat("left term", normfp, "\tright term\t", 
					normfk + con$btol * stepk * slopek, "descent", stepk * slopek, "\n")			

			#check Armijo condition
			if(normfp <= normfk + con$btol * stepk * slopek)
			{
				break
			}
		}
		
		inner.iter <- inner.iter + 1	
		stepk <- con$sigma * stepk
	}
	if(checkint(xk + stepk*dk, dimx, dimlam))
		normfp <- merit(xk + stepk*dk, dimx, dimlam, ...)
	else
		normfp <- Inf
	if(con$trace >= 3)
		cat("xk", xk, "\ndk", dk, "stepk", stepk, "normfk", normfk, "normfp", normfp, "\n")
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)

}			



# quadratic backtracking line search - minstep based on Dennis and Schnabel
linesearch.quad <- function(xk, dk, slopek, con, merit, ...)
{
	stepk <- 1			
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
	
	normfk <- merit(xk, ...)
	normfp <- merit(xk + stepk*dk, ...)
			
	while(stepk > minstep)
	{
		normfp <- merit(xk + stepk*dk, ...)
		
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")			
		
		
		#check Armijo condition
		if(normfp <= normfk + con$btol * stepk * slopek)
			break
		
		inner.iter <- inner.iter + 1	
		stepk <- - as.numeric( (stepk)^2 * slopek / 2 / (normfp - normfk - slopek*stepk)	)	

	}	
	if(con$trace >= 3)
		cat("xk", xk, "\ndk", dk, "stepk", stepk, "normfk", normfk, "normfp", normfp, "\n")
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		normfp=normfp, normfk=normfk)
}			
	
# quadratic backtracking line search customized for interior point methods
linesearch.quad.cond <- function(xk, dk, slopek, con, merit, checkint, dimx, dimlam, ...)
{
	stop("not yet implemented.")
	
	stepk <- 1
	inner.iter <- 0
	minstep <- con$xtol / max( abs(dk) / pmax(xk, 1) )		
	
	normfk <- merit(xk, dimx, dimlam, ...)	
	
	while(stepk > minstep)
	{
		#traces in R console	
		if(con$trace >= 3)
			cat("i", inner.iter, "\tstep", stepk, "\n")			
		
		if(checkint(xk + stepk*dk, dimx, dimlam))
		{			
			normfp <- merit(xk + stepk*dk, dimx, dimlam, ...)
		
			#traces in R console	
			if(con$trace >= 3)
				cat("i", inner.iter, "\tstep", stepk, "\tleft term", normfp, "\tright term\t", normfk + con$btol * stepk * slopek, "\n")
		
			#check Armijo condition
			if(normfp <= normfk + con$btol * stepk * slopek)
				break			
		}
		
		inner.iter <- inner.iter + 1	
		stepk <- - as.numeric( (stepk)^2 * slopek / 2 / (normfp - normfk - slopek*stepk) )	
		
	}	
	if(checkint(xk + stepk*dk, dimx, dimlam))
		normfp <- merit(xk + stepk*dk, dimx, dimlam, ...)
	else
		normfp <- Inf
	if(con$trace >= 3)
		cat("xk", xk, "\ndk", dk, "stepk", stepk, "normfk", normfk, "normfp", normfp, "\n")
	
	list(stepk=stepk, xk=xk, dk=dk, slopek=slopek, inner.iter=inner.iter, 
		 normfp=normfp, normfk=normfk)
}			

