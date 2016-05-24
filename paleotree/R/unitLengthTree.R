#' Scale Tree to Unit-Length
#' 
#' Rescales all edges of a phylogeny to be equal to 1 ("unit-length").
#' 
#' @details Probably not a good way to scale a tree for comparative studies.
#' What does it mean to scale every edge of the phylogeny to the same length?

#' @param tree an object of class phylo
#' @return Returns the modified phylogeny as an object of class phylo. Any
#' $root.time element is removed.
#' @seealso As an alternative to using unitLengthTree in comparative studies,
#' see \code{\link{timePaleoPhy}}
#' 
#' See also \code{speciationalTree} in the package geiger, which does
#' essentially the same thing as unitLengthTree.
#' @examples
#' 
#' set.seed(444)
#' tree <- rtree(10)
#' layout(1:2)
#' plot(tree)
#' plot(unitLengthTree(tree))
#' layout(1)
#' 
#' @export unitLengthTree
unitLengthTree<-function(tree){
	if(!inherits(tree, "phylo")){
		stop("tree is not of class phylo")
		}
	tree$edge.length<-rep(1,Nedge(tree))
	tree$root.time<-NULL
	return(tree)
	}
