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
### utility functions for semismooth Reformulation in GNE
###
###         R functions
### 


#functions of the SemiSmooth Reformulation of the GNEP
#z = (x , lambda)


#function phi of the SSR
funSSRcheck <- function(z, dimx, dimlam,
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	compl, argcompl, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	echo=FALSE)
{
	arg <- testargfunSSR(z, dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
						 compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint, echo)
	
	dimx <- arg$dimx
	dimlam <- arg$dimlam
	dimmu <- arg$dimmu	
	n <- sum(arg$dimx)
	m <- sum(arg$dimlam)
	nplayer <- length(arg$dimx)
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	mu <- z[-(1:(n+m))]

	
	#1st row is the begin index, 2nd row the end index
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) )

	GrLagri <- function(i) 
	{
		f <- function(j)
		{
			res <- arg$grobj(x, i, j, arg$arggrobj) 
			if(!is.null(arg$grconstr))
				res <- res + lam[index4lam[1,i]:index4lam[2,i]] %*% arg$grconstr(x, i, j, arg$arggrconstr)
			if(!is.null(arg$grjoint))
				res <- res + mu %*% arg$grjoint(x, j, arg$arggrjoint)
			res
		}
		#i index for player, j index for x_ij
		sapply(index4x[1,i]:index4x[2,i], f) 
	}
	complparti <- function(i)
	{
		#i index for player
		Constri <- arg$constr(x, i, arg$argconstr)
		arg$compl(-Constri, lam[index4lam[1,i]:index4lam[2,i]], arg$argcompl)
	}
	
	
	part1 <- sapply(1:nplayer, GrLagri)
	if(!is.null(arg$grconstr))
		part2a <- sapply(1:nplayer, complparti)
	else
		part2a <- NULL
	if(!is.null(arg$joint))
		part2b <- arg$compl(-arg$joint(x, arg$argjoint), mu, arg$argcompl)
	else
		part2b <- NULL
	
	unlist( c(part1, part2a, part2b) )
}









#Jacobian of phi of the SSR
jacSSRcheck <- function(z, dimx, dimlam, 
	heobj, argheobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	gcompla, gcomplb, argcompl, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	hejoint, arghejoint,
	echo=FALSE)
{
	
	arg <- testargjacSSR(z, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
		heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint,
		hejoint, arghejoint, echo=FALSE)
	
	dimx <- arg$dimx
	dimlam <- arg$dimlam
	dimmu <- arg$dimmu
	n <- sum(dimx)
	m <- sum(dimlam)
	nplayer <- length(dimx)
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	mu <- z[-(1:(n+m))]
	
	#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) )
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )	
	
	grjgriLagri <- function(i, j)
	{
		f <- function(k)
		{
			res <- arg$heobj(z, i, j, k, arg$argheobj) 
			if(!is.null(arg$heconstr))
				res <- res + lam[index4lam[1,i]:index4lam[2,i]] %*% arg$heconstr(z, i, j, k, arg$argheconstr)
			if(!is.null(arg$hejoint))
				res <- res + mu %*% arg$hejoint(x, j, k, arg$arghejoint)
			res
		}
		
		#i index for player, j index for x_j, k for x_i_k		
		sapply(index4x[1,i]:index4x[2,i], f)
	}
	grjConstri <- function(i) 
	{
		#i index for player
		sapply(index4x[1,i]:index4x[2,i], function(j) arg$grconstr(z, i, j, arg$arggrconstr))
	}
	grjJointi <- function(i) 
	{
		#i index for player
		sapply(index4x[1,i]:index4x[2,i], function(j) arg$grjoint(z, j, arg$arggrjoint))
	}
	
	grxjcomplparti <- function(i, j)
	{
		#i index for player, j index for x_j		
		GrjConstri <- arg$grconstr(z, i, j, arg$arggrconstr)
		Constri <- arg$constr(z, i, arg$argconstr)
		- GrjConstri * arg$gcompla(-Constri, lam[index4lam[1,i]:index4lam[2,i]], arg$argcompl)
	}
	grlamicomplparti <- function(i)
	{
		#i index for player
		Constri <- arg$constr(z, i, arg$argconstr)
		arg$gcomplb(-Constri, lam[index4lam[1,i]:index4lam[2,i]], arg$argcompl)
	}
	
	grxjcomplpartb <- function(j)
	{
		#i index for player, j index for x_j		
		GrjConstr <- arg$grjoint(z, j, arg$arggrjoint)
		Constr <- arg$joint(z, arg$argjoint)
		- GrjConstr * arg$gcompla(-Constr, mu, arg$argcompl)
	}	
	grmucomplpartb <- function()
	{
		#i index for player
		Constr <- arg$joint(z, arg$argjoint)
		arg$gcomplb(-Constr, mu, arg$argcompl)
	}	
	
	#Hessian matrix of the Lagrangian
	ggL <- matrix(0, sum(dimx), sum(dimx))
	for(i in 1:nplayer)
		ggL[index4x[1,i]:index4x[2,i] , ] <- sapply(1:sum(dimx), function(j) grjgriLagri(i,j))	
	#gradient of the constraint function
	gG <- matrix(0, sum(dimx), sum(dimlam))
	if(!is.null(arg$heconstr))
	for(i in 1:nplayer)
		gG[index4x[1,i]:index4x[2,i] , index4lam[1,i]:index4lam[2,i]] <- t( grjConstri(i) )
	#gradient of the constraint function
	gH <- matrix(0, sum(dimx), dimmu)
	if(!is.null(arg$hejoint))
	for(i in 1:nplayer)
		gH[index4x[1,i]:index4x[2,i] , ] <- t( grjJointi(i) )
	
	
	#gradient of the complementarity function
	gxcompla <- matrix(0, sum(dimlam), sum(dimx))
	if(!is.null(arg$heconstr))
	for(i in 1:nplayer)
	for(j in 1:sum(dimx))
		gxcompla[index4lam[1,i]:index4lam[2,i] , j] <-  grxjcomplparti(i,j)
	
	#gradient of the complementarity function
	glcompla <- matrix(0, sum(dimlam), sum(dimlam))
	if(!is.null(arg$heconstr))
	for(i in 1:nplayer)
		glcompla[index4lam[1,i]:index4lam[2,i] , index4lam[1,i]:index4lam[2,i]] <- diag(grlamicomplparti(i), dimlam[i], dimlam[i])
	
	m01 <- matrix(0, sum(dimlam), dimmu)
	
	#gradient of the complementarity function
	gxcomplb <- matrix(0, dimmu, sum(dimx))
	if(!is.null(arg$hejoint))
	for(j in 1:sum(dimx))
		gxcomplb[, j] <-  grxjcomplpartb(j)
	
	#gradient of the complementarity function
	if(!is.null(arg$hejoint))
	glcomplb <- diag(grmucomplpartb(), dimmu, dimmu)
	
	m02 <- matrix(0, dimmu, sum(dimlam))
	
	if(!is.null(arg$hejoint) && !is.null(arg$heconstr))
		return(rbind(cbind(ggL,			gG,			gH),
					 cbind(gxcompla,	glcompla,	m01),
					 cbind(gxcomplb,	m02,		glcomplb)) )
	else if(is.null(arg$hejoint) && !is.null(arg$heconstr))
		return(rbind(cbind(ggL,		gG),
					 cbind(gxcompla,glcompla)) )
	else if(!is.null(arg$hejoint) && is.null(arg$heconstr))
		return(rbind(cbind(ggL,			gH),
					 cbind(gxcomplb,	glcomplb)) )
	else
		return(ggL)
	
}




