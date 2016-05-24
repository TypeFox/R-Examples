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


copula <- function(expdata) # symetric MI matrix is input
{
numofcol <- ncol(expdata)
numofrow <- nrow(expdata)
a <- matrix(c(0),numofrow,numofcol)
aa <- matrix(c(0),numofrow,numofcol)
a <- expdata
 for(i in 1:numofrow) {
  aa[i,] <- matrix( rank(a[i,], ties.method= "first"),ncol=numofcol )
  			}       
b<-(aa- 0.5)/numofcol 
rownames(b) <- rownames(expdata)
colnames(b) <- colnames(expdata)     

b
}

