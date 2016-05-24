#'Calculates maxdist for all pairs of taxa in a phylogeny.   
#'
#'allmaxdist Uses maxdist to calcualte the maximum phenotypic.  distance between the ancestors of all pairs of taxa in a phylogeny.  By default outputs these as a matrix, but can also output a list.  Can take some time to run for large trees.    
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param mat Whether or not to export the values in a matrix (default) or a list
#'
#'@details  Regarding the output:  the matrix is better organized, but the list avoids all the zeroes and is probably better for making histograms.
#'
#'@return A matrix or list of all maxdist values for all pairs of taxa in the phylogeny.  
#'
#'@import ape geiger MASS phytools
#'
#'@export
#'
#'@references Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
#'and evolution in R langauge. Bioinformatics, 20, 289-290.
#'
#'Revell, L. J. (2012) phytools: An R package for phylogenetic comparative 
#'biology (and other things). Methods Ecol. Evol. 3 217-223.
#'
#'@examples
#'
#'phyl<-rtree(10)
#'phendata<-fastBM(phyl,nsim=2)
#'answer<-allmaxdist(phyl,phendata,mat=TRUE)


allmaxdist<-function(phyl,phendata,mat=TRUE)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    

#The function

dims<-dim(phyl$edge)

ntips<-(dims[1]/2)+1

allmaxes<-matrix(0,ntips,ntips)

maxeslist<-c()

t1<-1

while (t1<ntips) {

	t2<-t1+1

	while (t2<=ntips) {

		m<-maxdist(phyl,phendata,t1,t2)

		allmaxes[t1,t2]<-m

		maxeslist<-c(maxeslist,m)

		t2<-t2+1

		} 

	t1<-t1+1

	}

if (mat==TRUE) {answer<-allmaxes} else{answer<-maxeslist}

}

