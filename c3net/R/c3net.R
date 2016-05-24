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


c3net <- function(dataset, cop=TRUE, alpha=0.01, methodstep1="cutoff", cutoffMI= 0, MTCmethod="BH", itnum=5, network=FALSE)
{
      net <- NULL
      if(cop==TRUE){
      dataset <- copula(dataset)
      }
      mim<-makemim(dataset)
      
      if( methodstep1=="cutoff" ) {
	 if(cutoffMI== 0) Ic <- mean(mim[upper.tri(mim)])
	  else Ic <- cutoffMI	

	mim[mim < Ic] <-0 
	net <- c3(mim)
	}
      else if( methodstep1=="MTC") {
	res <- sigtestMTC(dataset, alpha, itnum, methodsig=MTCmethod)
        net <- c3(res$Inew) 
	}	 
      else {
	res <- sigtestp(dataset, alpha, itnum)
        net <- c3(res$Inew) 
	}

if(network==TRUE) netplot(net)
 
net 
}
