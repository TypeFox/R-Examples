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


sigtestMTC <- function(data, alpha, itnum, methodsig="BH")   
{



################# 
res <- sigtestp(data, alpha, itnum) 	
################# 
mim <- res$mim
### FDR process starts

Inew <- mim
ngene<- ncol(mim)

MIs <- mim[upper.tri(mim, diag = FALSE)]

vg <- res$vg

pMI <- matrix(c(1),ngene,ngene)
sum_vg <- sum(vg)
vgnum <- length(vg)
E <- length(MIs)

pvg <- c()
max_vg <- max(vg)

for(i in 1: (ngene-1)) for(ix in (i+1):ngene)
{
j <- 0
ii <- 1

if(max_vg <= mim[i,ix])
{
 pMI[i,ix] <- 0 
 pMI[ix,i] <- 0 
}

if(max_vg > mim[i,ix])
{

ind <- which( vg > mim[i,ix])
LL <- length(vg)
ptemp <- length(ind)/LL
 pMI[i,ix] <- ptemp
 pMI[ix,i] <- ptemp
} 

}

p_val <- pMI[upper.tri(pMI, diag = FALSE)]

	padj <- p.adjust(p_val,  method = methodsig)
	MIedges <- MIs 
	for (jj in 1:E) { if(padj[jj] > alpha) MIs[jj] <- 0 }  #eliminate nonsignificant MIs 

	Inew[upper.tri(Inew)] <- MIs
	txI<-t(Inew)
	Inew[lower.tri(Inew)] <- txI[lower.tri(txI)]
	diag(Inew) <- 0


pMIadj <- matrix(c(1),ngene,ngene)
pMIadj[upper.tri(pMIadj)] <- padj


	res <- new.env()  
	assign("I0", res$I0, envir=res)
	assign("vg", res$vg, envir=res)
	assign("padj", padj, envir=res)
	assign("p_val", p_val, envir=res)
	assign("pMIadj", padj, envir=res)
	assign("pMI", p_val, envir=res)
	assign("Inew", Inew, envir=res)
	assign("mim", mim, envir=res)
	assign("MIedges", MIedges, envir=res)

res   

} 


