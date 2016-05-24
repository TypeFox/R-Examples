#'Calculates the maximum phenotypic distance between the lineages leading to a pair of taxa. 
#'
#'maxdist uses ancestral state reconstruction to determine the maximum distance between any ancestors of those two taxa.   
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param phendata Phenotypic data for all tips
#'@param t1 The first taxon of interest
#'@param t2 The second taxon of interest
#'
#'@details Returns the maximum Euclidean distance between any pair of ancestors of the two taxa, whether or not those two ancestors are contemporaries.  
#'
#'@return The maximum phenotypic distance between the two taxa
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
#'answer<-maxdist(phyl,phendata,1,10)


maxdist<-function(phyl,phendata,t1,t2)

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

if (is.finite(t1)) {}
	else {t1<-labelstonumbers(phyl,t1)}

if (is.finite(t2)) {}
	else {t2<-labelstonumbers(phyl,t2)}

if (t1>length(phyl$tip))
	stop("your first tip isn't in the phylogeny.")

if (t2>length(phyl$tip))
	stop("your second tip isn't in the phylogeny.")




#The function

alldata<-multianc(phyl,phendata)

anctimes<-node.depth.edgelength(phyl)

combineddata<-cbind(anctimes,alldata)

#Then we go through and grab just the mrca of our two tips, and all nodes 
#between the mrca and the tips.  This part will be based on findanc.

mrcas<-mrca(phyl)

mrcat1t2<-mrcas[t1,t2]

#First, the part for t1

t1path<-combineddata[t1 ,]

anc<-findanc(phyl,t1)

t1path<-rbind(t1path, combineddata[anc[1] ,])

while (anc[1] != mrcat1t2) {

	anc<-findanc(phyl, anc[1])

	t1path<-rbind(t1path, combineddata[anc[1] ,])

	}

#And then something identical for t2

t2path<-combineddata[t2 ,]

anc<-findanc(phyl,t2)

t2path<-rbind(t2path, combineddata[anc[1] ,])

while (anc[1] != mrcat1t2) {

	anc<-findanc(phyl, anc[1])

	t2path<-rbind(t2path, combineddata[anc[1] ,])

	}

#We just want the trait values for now

dims<-dim(t1path)

t1pathtraits<-t1path[, 2:dims[2]]
t2pathtraits<-t2path[, 2:dims[2]]

#And we'll bind them together

alltraits<-rbind(t1pathtraits, t2pathtraits)

phendist<-dist(alltraits,method="euclidean",diag=TRUE, upper=TRUE)

answer<-max(phendist)

}
