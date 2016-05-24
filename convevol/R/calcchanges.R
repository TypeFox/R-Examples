#'Calculates all phenotypic changes that occur on all branches of a phylogeny.  
#'
#'Calculates the Euclidean distance between all ancestors and descendants on a phylogeny to reconstruct the phenotypic changes that occur along all edges of a phylogeny.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata A matrix of phenotypic data for all tips, with taxa in rows and characters in columns.
#'
#'@details Calculates the Euclidean distance between all ancestors and descendants on a phylogeny to reconstruct the phenotypic changes that occur along all edges of a phylogeny. 
#'
#'@return A vector in which each element represents an edge of the phylogeny, and the values are the magnitudes of evolutionary change that occur along those edges.
#'
#'import ape geiger phytools
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
#'phendata<-fastBM(phyl,nsim=5)
#'changes<-calcchanges(phyl,phendata)

calcchanges<-function(phyl,phendata)

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (nrow(phendata) != length(phyl$tip)) 
	stop("your data matrix must have the same number of rows as tips in the tree.")
    
if (is.null(rownames(phendata))) {
	warning("no row names for data.  Assuming that the rows are in the same order as tips.")
      rownames(X) <- phyl$tip.label
	}

#The function.  First, ancestral states are reconstructed for the phenotypic data.  


allvals<-multianc(phyl,phendata)  

pdims<-dim(phyl$edge)

changes<-rep(0,pdims[1])

j<-1

for (j in 1:pdims[1]) {

	ancnode<-phyl$edge[j,1]
	desnode<-phyl$edge[j,2]

	ancvals<-allvals[ancnode ,]
	desvals<-allvals[desnode ,]

	change<-(sum((ancvals-desvals)^2))^0.5

	changes[j]<-change

	}

answer<-changes

}