#This file belongs to
#c3net: C3NET, <https://r-forge.r-project.org/projects/c3net/>
#This R package allows inferring regulatory networks from expression data using C3NET.
#The inferred network consists of only direct physical interactions.
## Copyright (C) January 2011 Gokmen Altay <altayscience@gmail.com>
## This program is a free software for only academic useage but not for commercial useage; you can redistribute it and/or
## modify it under the terms of the GNU GENERAL PUBLIC LICENSE
## either version 3 of the License, or any later version.
##
## This program is distributed WITHOUT ANY WARRANTY; 
## You can get a copy of the GNU GENERAL PUBLIC LICENSE
## from
## http://www.gnu.org/licenses/gpl.html
## See the licence information for the dependent package from
## igraph package itself.


makemim <- function(expdata) # symetric MI matrix is input
{
gnames<-rownames(expdata)
txdata <-t(expdata)
genenum <- ncol(txdata)
	rho2 <- cor(txdata)^2
	for(irh in 1:genenum) for(jrh in 1:genenum) 
		{ if(rho2[irh,jrh]>=1) {rho2[irh,jrh] <- 0.9999999} 
		}
	mim <- -0.5*log(1-rho2)
	mim <- abs(mim)  
	diag(mim)<-0
	rownames(mim) <-gnames
	colnames(mim) <-gnames

mim
}


