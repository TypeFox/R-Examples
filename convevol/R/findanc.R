#'Find the ancestor of a given node in a phylogeny
#'
#'This function will find the ancestor of a given node in a phylogeny.  It will return a two-element vector, which will contain both the node of the ancestor and the number of the edge that connects the node and ancestor.
#'
#'@param phyl The phylogeny of interest in phylo format
#'@param node The number of the node that you want the ancestor for
#'
#'@details   Returns a two-element vector.  The node of the ancestor is first;  the edge that connects that node with its ancestor is second.  
#'
#'@return A two-element vector, where the first element is the node of the ancestor and the second element is the number of the edge that connects the node and ancestor (i.e., the row number in phyl$edge).
#'
#'@import ape geiger MASS phytools
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
#'ancestor<-findanc(phyl,1)



findanc<-function(phyl,node)

#This function will find the ancestor of a given node in a phylogeny.  And 
#that's all it'll do.  phyl is the phylogeny and node is the number of the 
#node that we're interested in.

#It will return a two-element vector, which will contain both the node of the 
#ancestor, and the number of the branch that connects the node and ancestor.

{

#Error checking

if (class(phyl) != "phylo") 
	stop("your tree must be class 'phylo.'")

if (node>dim(phyl$edge)[1]+1)
	stop("there are not that many nodes in the tree")

#The function

if(is.character(node)==TRUE){node<-labelstonumbers(phyl,node)}

i<-1

tempnode<-phyl$edge[i,2]

while (tempnode != node) {

	i<-i+1

	tempnode<-phyl$edge[i,2]
	tempnode
	}
        
anc<-phyl$edge[i,1]

answer<-c(anc,i)
}

