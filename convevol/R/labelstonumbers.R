#'Converts taxon names to tip/edge numbers
#'
#'Converts taxon names to corresponding tip/edge numbers in the phylogeny.  
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param tips The names of the tips in question
#'
#'@details Simply reads in taxon names, determines which tip or edge number they correspond to, and returns those values   
#'
#'@return The numbers of all of the tips of interest.  
#'
#'@import ape
#'
#'@export
#'
#'@references Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics
#'and evolution in R langauge. Bioinformatics, 20, 289-290.
#'Paradis, E. (2012) Analysis of Phylogenetics and Evolution with R (Second Edition). New York: Springer. 
#'
#'@examples
#'
#'phyl<-rtree(10)
#'nums<-labelstonumbers(phyl,c("t1","t2","t3"))


labelstonumbers<-function(phyl,tips)

#This will take a set of tip names, find the corresponding tip number (for use 
#in phyl$edge, for example), and return them in order.  

{

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

newtips<-c()


for(j in 1:length(tips)) {

	for(i in 1:length(phyl$tip.label)) {
	
		if(phyl$tip.label[i]==tips[j]) {newtips<-c(newtips,i)}

		}

	}

newtips

}

