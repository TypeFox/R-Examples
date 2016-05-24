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
#z = (x, lam, w)
#with size (n, m, m)


funCER <- function(z, dimx, dimlam, 
	grobj, arggrobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	echo=FALSE)
{
	arg <- testargfunCER(z, dimx, dimlam, grobj, arggrobj, constr, argconstr, grconstr, arggrconstr, 
						 dimmu, joint, argjoint, grjoint, arggrjoint, echo)
	
	dimx <- arg$dimx
	dimlam <- arg$dimlam
	dimmu <- arg$dimmu	
	n <- sum(arg$dimx)
	m <- sum(arg$dimlam)
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	mu <- z[(n+m+1):(n+m+dimmu)]
	w <- z[-(1:(n+m+dimmu))]
	
	part1 <- funSSR(z[1:(n+m+dimmu)], dimx, dimlam, grobj, arggrobj, constr, argconstr, grconstr, arggrconstr, 
					compl=phiFB, argcompl=NULL, dimmu, joint, argjoint, grjoint, arggrjoint)[1:n]

	Constri <- function(i) arg$constr(z, i, arg$argconstr)
	if(!is.null(arg$grconstr))	
	{
		part2a <- unlist(sapply(1:arg$nplayer, Constri)) + w[1:m]
		part3a <- lam * w[1:m]
	}else
		part2a <- part3a <- NULL
		
	if(!is.null(arg$joint))	
	{
		part2b <- arg$joint(x, arg$argjoint) + w[m+1:dimmu]
		part3b <- mu * w[m+1:dimmu]
	}else
		part2b <- part3b <- NULL

	c( part1, part2a, part2b, part3a, part3b )
}



#z = (x, lam, w)
#with size (n, m, m)
jacCER <- function(z, dimx, dimlam,
	heobj, argheobj, 
	constr, argconstr,  
	grconstr, arggrconstr, 
	heconstr, argheconstr,
	dimmu, joint, argjoint,
	grjoint, arggrjoint,
	hejoint, arghejoint,
	echo=FALSE)
{
	arg <- testargjacCER(z, dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
						heconstr, argheconstr, dimmu, joint, argjoint, grjoint, arggrjoint, 
						 hejoint, arghejoint, echo)
									 
	
	dimx <- arg$dimx
	dimlam <- arg$dimlam
	dimmu <- arg$dimmu	
	n <- sum(arg$dimx)
	m <- sum(arg$dimlam)
	p <- dimmu
	x <- z[1:n]
	lam <- z[(n+1):(n+m)]
	mu <- z[(n+m+1):(n+m+p)]
	w1 <- z[(n+m+p+1):(n+m+p+m)]
	w2 <- z[(n+2*m+p+1):(n+2*m+2*p)]
	nplayer <- arg$nplayer
	

#1st row is the begin index, 2nd row the end index
	index4lam <- rbind( cumsum(dimlam) - dimlam + 1, cumsum(dimlam) )
	index4x <- rbind( cumsum(dimx) - dimx + 1, cumsum(dimx) )	
	
	GrjConstri <- function(i) 
	{
		#i index for player
		sapply(index4x[1,i]:index4x[2,i], function(j) arg$grconstr(z, i, j, arg$arggrconstr))
	}
	jacconstrij <- function(i, j)
	{
		#i index for player, j index for x_j		
		arg$grconstr(z, i, j, arg$arggrconstr)
	}
	jacjointj <- function(j)
	{
		#j index for x_j		
		arg$grjoint(z, j, arg$arggrjoint)
	}
	
	partSSR <- jacSSR(z[1:(n+m+p)], dimx, dimlam, heobj, argheobj, constr, argconstr, grconstr, arggrconstr, 
					  heconstr, argheconstr, gcompla=GrAphiFB, gcomplb=GrBphiFB, argcompl=NULL, 
					  dimmu, joint, argjoint, grjoint, arggrjoint, hejoint, arghejoint)
	
#Hessian matrix of the Lagrangian
	ggL <- partSSR[1:n, 1:n]

	
#gradient of the constraint function
	if(!is.null(arg$heconstr))
		gG <- partSSR[1:n, n+ 1:m] 
		
#gradient of the joint function
	if(!is.null(arg$hejoint))	
		gH <- partSSR[1:n, n+m+ 1:p]
	
#Jacobian of the constraint function
	jacG <- matrix(0, m, n)
	if(!is.null(arg$heconstr))	
	for(i in 1:nplayer)
	for(j in 1:sum(dimx))
		jacG[index4lam[1,i]:index4lam[2,i] , j] <- jacconstrij(i,j)
	
#Jacobian of the joint function	
	jacH <- matrix(0, p, n)
	if(!is.null(arg$hejoint))	
	for(j in 1:sum(dimx))
		jacH[, j] <- jacjointj(j) 	
	
	m0mm <- matrix(0, m, m)
	m0pp <- matrix(0, p, p)
	m0mp <- matrix(0, m, p)
	m0mn <- matrix(0, m, n)
	m0pn <- matrix(0, p, n)
	
	if(!is.null(arg$hejoint) && !is.null(arg$heconstr))
		res <- rbind(cbind(ggL,	gG,			gH,		t(m0mn),	t(m0pn)),
					cbind(jacG,	m0mm,		m0mp,	diag(m),	m0mp),
					cbind(jacH,	t(m0mp),	m0pp,	t(m0mp),	diag(p)),
					cbind(m0mn, diag(w1),		m0mp,	diag(lam),	m0mp),
					cbind(m0pn, t(m0mp),	diag(w2),	t(m0mp),	diag(mu))
					 )
	else if(is.null(arg$hejoint) && !is.null(arg$heconstr))
		res <- rbind(cbind(ggL,	gG,			t(m0mn)),
				 cbind(jacG,	m0mm,		diag(m)),
				 cbind(m0mn,	diag(w1),	diag(lam))
				 )	
	else if(!is.null(arg$hejoint) && is.null(arg$heconstr))
		res <- rbind(cbind(ggL,	gH,			t(m0pn)),
				 cbind(jacH,	m0pp,		diag(p)),
				 cbind(m0pn,	diag(w2),	diag(mu))
				 )	
	else
		res <- ggL
	return(res)
}


