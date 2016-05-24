#' Split Tip Taxa by Root Divergence
#' 
#' Sorts terminal taxa into groups descended from each lineage splitting off of
#' the root node.
#' 
#' @details This function can be useful for studying the timing in the order
#' of appearance of descended from different lineages descended from the first
#' bifurcation.
#' 
#' @param tree A phylo object
#' @return Returns a list with each element a character vector containing the
#' names of terminal taxa descended from each lineage splitting off of the root
#' node.
#' @author David W. Bapst
#' @examples
#' 
#' tree<-rtree(100)
#' rootSplit(tree)
#' 
#' @export rootSplit
rootSplit<-function(tree){
	#returns a list with the daughter taxa of the two clades at the root split
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	tips<-lapply(tree$edge[tree$edge[,1]==(Ntip(tree)+1),2],function(zz) 
		if(zz>Ntip(tree)){tree$tip.label[prop.part(tree)[[zz-Ntip(tree)]]]
			}else{tree$tip.label[zz]})
	return(tips)
	}
