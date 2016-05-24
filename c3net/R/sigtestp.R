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


sigtestp <- function(data, alpha, itnum)   
{
	gnames<-rownames(data)
	a<-as.matrix(data)
	numbrow <- nrow(a)
	numbcol <- ncol(a)
	
	numbdataset <- itnum   #number of exp random matrices	
	v<-c()    #final MI vector 

	for(i in 1:numbdataset)
	{
	expshufled<-matrix(sample(a), numbrow,numbcol) ##randomize all elements
	
	######################## copula transform
	expdata<-copula(expshufled)  
	
	########################
	mim <- makemim(expdata)
	###################

	c<-mim
	ctry<-c
	ccc<-ctry[upper.tri(ctry)]
	ccc<-ccc[ccc!=0]

	t1<-as.vector(ccc)
	v<-c(v,t1)
	} #end for i
	vg<-sort(v)

	
################# Sample distribution and threshold function  ends ##########################

### obtain MIM matrix
	######################## copula transform
	expdata<-copula(data) 
	
	########################
	mim <- makemim(expdata)
	###################
	rownames(mim) <-gnames
	colnames(mim) <-gnames

######


## from alpha, MI threshold estimation
L <- length(vg)
temp <- ceiling((1-alpha)*L)
I0 <- vg[temp]
#############
Inew <- mim
Inew[Inew < I0] <-0 

	res <- new.env()  
	assign("I0", I0, envir=res)
	assign("vg", vg, envir=res)
	assign("Inew", Inew, envir=res)
	assign("mim", mim, envir=res)


res   

} 


