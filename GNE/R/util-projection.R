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
### utility functions for projection in GNE
###
###         R functions
### 




projector <- function(z, g, jacg, bounds=c(0, 10), echo=FALSE, ...)
{
	merit <- function(y, z) 1/2 * sum( (y - z)^2 )
	grmerit <- function(y, z) y - z	
	
	zinit <- rejection(g, length(z), bounds[1], bounds[2], ...)
	
	
	if(echo)
		cat("init", zinit,"\nz", z, "\n")
	

#	print(merit(zinit, z))
#	print(grmerit(zinit, z))
#	print(-g(zinit, ...))
#	print(-jacg(zinit, ...))
	
	res <- constrOptim.nl(zinit, 
				fn=function(y, ...) merit(y, z), 
				gr=function(y, ...) grmerit(y, z), 
				hin=function(y, ...) -g(y, ...), 
				hin.jac=function(y, ...) -jacg(y, ...), 
				control.outer=list(trace=FALSE),
				control.optim=list(trace=FALSE), ...)
	
	if(echo)
		print(res)
	
	res$par
}

