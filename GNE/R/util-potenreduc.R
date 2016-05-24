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
### utility functions for potential reduction algorithm in GNE
###
###         R functions
### 



#functions of the potential reduction algorithm



#logarithm potential function	
potential.ce <- function(u, n, zeta)
{
	res <- zeta * log( sum(u^2) )
	if(length(u) > n)
	{
		if(any(u[-(1:n)] < 0))
			stop(paste("u has negative components:", paste(u, collapse=", "), "\n"))
		res <- res - sum( log( u[-(1:n)] ) )
	}
	res
}

#gradient of logarithm potential function		
gradpotential.ce <- function(u, n, zeta)	
{
	normu <- sum(u^2)
	if(length(u) > n)
	{
		return( c( 2*zeta / normu * u[1:n], 2*zeta / normu * u[-(1:n)] - 1/u[-(1:n)] ) )
	}else
	{
		return( 2*zeta / normu * u ) 
	}
}	



#merit function for the constrained equation	
#z = (x, lam, w)
#with size (n, m, m)
psi.ce <- function(z, dimx, dimlam, Hfinal, argfun, zeta)
	as.numeric( potential.ce( Hfinal(z, argfun=argfun), sum(dimx), zeta) )

#gradient of the merit function	
#z = (x, lam, w)
#with size (n, m, m)
gradpsi.ce <- function(z, dimx, dimlam, Hfinal, jacHfinal, 
argfun, argjac, zeta)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	w <- z[(n+m+1):(n+2*m)]
	
	Hz <- Hfinal(z, argfun=argfun)
	Jz <- jacHfinal(z, argjac=argjac)
	
	
#t(Jz) =
#	t(Jz[1:n, 1:n])				t(Jz[(n+1):(n+m), 1:n])		0
#	t(Jz[1:n, (n+1):(n+m)])		0							diag(w)
#	0							Identity_m					diag(lam)			
	
	grpHz <- gradpotential.ce(Hz, n, zeta)
	
	if(m>0)
	{
		gradpsi <- c(	t(Jz[1:n, 1:n]) %*% grpHz[1:n] + t(Jz[(n+1):(n+m), 1:n]) %*% grpHz[(n+1):(n+m)], 
				 t(Jz[1:n, (n+1):(n+m)])  %*% grpHz[1:n] + diag(w) %*% grpHz[(n+m+1):(n+2*m)], 
				 grpHz[1:m+n] + diag(lam) %*% grpHz[(n+m+1):(n+2*m)])
	}else
	{
		gradpsi <- t(Jz[1:n, 1:n]) %*% grpHz[1:n]
	}
	return(gradpsi)
}


#check interior
#z = (x, lam, w)
#with size (n, m, m)
checkint.ce <- function(z, dimx, dimlam)
{
	n <- sum(dimx)
	m <- sum(dimlam)
	
	if(m > 0)
	{
		lam <- z[(n+1):(n+m)]
		w <- z[(n+m+1):(n+2*m)]
		return( all(lam > 0, w > 0) )
	}else
		return(TRUE)
}

