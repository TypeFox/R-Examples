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
### utility functions for trust-region in GNE
###
###         R functions
### 


# trust region utility functions

#Newton point : pnewton of 2-norm npnewton
#Cauchy point (steepest descent direction) : pcauchy of 2-norm npcauchy
#trust region radius : del
poweldogleg <- function(n, pnewton, npnewton, pcauchy, npcauchy, delta)
{
	if(length(pcauchy) != n || length(pnewton) != n)
		stop("wrong argument in pdleg.")
	if(length(npcauchy) > 1 || length(npnewton) > 1)
		stop("wrong argument in pdleg.")
	
	if(npnewton <= delta) #Newton step smaller than trust region radius
	{
		nextp <- pnewton
		nexttype <- "Newton step"
	}else if(npcauchy >= delta) #Cauchy step greater than trust region radius
	{
		nextp <- delta / npcauchy * pcauchy
		nexttype <- "Cauchy step"
	}else #convex combination between Newton and Cauchy points
	{
		pc_pn_pc <- crossprod(pcauchy, pnewton-pcauchy)
		pn_pc2 <- crossprod(pnewton-pcauchy, pnewton-pcauchy)
		lam <- - pc_pn_pc + sqrt(pc_pn_pc^2  - pn_pc2*(npcauchy - delta))
		lam <- as.numeric(lam / pn_pc2)
		nextp <- lam * pnewton + (1-lam) * pcauchy
		nexttype <- "convex point"
	}
	
	list(nextp=nextp, type=nexttype)
}

#Newton point : pnewton of 2-norm npnewton
#Cauchy point (steepest descent direction) : pcauchy of 2-norm npcauchy
#lower / upper bounds : l, u
#scaling vector : scalvec
#control list : control
truncpoweldogleg <- function(n, pnewton, npnewton, pcauchy, npcauchy, 
	delta, l, u, scalvec, xk, control)
{
	if(length(pcauchy) != n || length(pnewton) != n)
		stop("wrong argument in pdleg.")
	if(length(npcauchy) > 1 || length(npnewton) > 1)
		stop("wrong argument in pdleg.")
	
	#scaled vectors
	scpnew <- scalvec * pnewton
	scpcau <- scalvec * pcauchy
	#scaled norm
	scnpnew <- sqrt(sum(scpnew^2))
	scnpcau <- sqrt(sum(scpcau^2))
	
	if(scnpnew <= delta) #Newton step smaller than trust region radius
	{
		nextp <- pnewton
		nexttype <- "Newton step"
		
		stepsize <- min(pmax((l - xk)/pnewton, (u - xk)/pnewton))
		if(stepsize <= 1)
		{	
			nextp <- max(control$theta, 1-npnewton) * stepsize * pnewton
			nexttype <- "truncated Newton step"
		}
	}else if(scnpcau >= delta) #Cauchy step greater than trust region radius
	{
		nextp <- delta / npcauchy * pcauchy
		nexttype <- "Cauchy step"

		stepsize <- min(pmax((l - xk)/pcauchy, (u - xk)/pcauchy))
		if(stepsize <= 1)
		{	
			nextp <- max(control$theta, 1-npcauchy) * stepsize * pcauchy
			nexttype <- "truncated Cauchy step"
		}
	}else #convex combination between Newton and Cauchy points
	{
		pc_pn_pc <- crossprod(scpcau, scnpnew-scpcau)
		pn_pc2 <- crossprod(scnpnew-scpcau, scnpnew-scpcau)
		lam <- - pc_pn_pc + sqrt(pc_pn_pc^2  - pn_pc2*(scnpcau - delta))
		lam <- as.numeric(lam / pn_pc2)
		nextp <-  (lam * scnpnew + (1-lam) * scpcau) / scalvec
		nnextp <- sqrt(sum(nextp^2))
		nexttype <- "convex point"
		
		stepsize <- min(pmax((l - xk)/nextp, (u - xk)/nextp))
		if(stepsize <= 1)
		{	
			nextp <- max(control$theta, 1-nnextp) * stepsize * nextp
			nexttype <- "truncated convex step"
		}
		
	}
	
	list(nextp=nextp, type=nexttype)
}


