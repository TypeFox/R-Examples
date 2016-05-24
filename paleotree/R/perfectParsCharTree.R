#' Simulate a Set of Parsimony-Informative Characters for a Phylogeny
#'
#' Creates a simulated set of parsimony-informative characters for a given rooted phylogeny,
#' with characters shared out equally across nodes in the phylogeny, with any remaining characters
#' assigned randomly to nodes.

#' @details
#' This function takes some a tree and places a number of binary characters on the tree, with character states arranged  
#' as if the derived condition was gained once, at a single node, and never lost. This ensures that the resulting simulated
#' character matrices have no character conflict, supporting a single solution under maximum parsimony. 
#'
#' If \code{nchar} is greater than the number of nodes on the input phylogeny (ignoring the root), then characters are first
#' placed to evenly cover all nodes, with as many full passes of tree as possible. Any characters in excess are placed at random
#' nodes, without replacement. In other words, if a tree has 10 nodes (plus the root) and 25 characters are simulated, 20 of those
#' characters will consist of two 10-character 'full passes' of the tree. The remaining five will be randomly dropped on the tree.
#'
#' If few characters are simulated than the number of nodes, these are randomly placed on the given topology without replacement,
#' just as described above. 
#'
#' This function assumes, like almost every function in paleotree, that the tree given is rooted, even if the
#' most basal node is a polytomy.

#' @param tree A phylogeny of class 'phylo'

#' @param nchar Number of parsimonious binary characters to simulate on the phylogeny.

#' @return A matrix of \code{nchar} parsimonious binary characters for each taxon on \code{tree}, with states 0 and 1.

#' @author David W. Bapst 

#' @examples
#' data(retiolitinae)
#'
#' #fewer characters than nodes
#' perfectParsCharTree(retioTree,nchar=10)
#'
#' #same as number of nodes (minus root)
#' perfectParsCharTree(retioTree,nchar=12)
#'
#' #more characters than the number of nodes
#' perfectParsCharTree(retioTree,nchar=20)

#' @export perfectParsCharTree
perfectParsCharTree<-function(tree,nchar){
	#checks
	if(!inherits(tree,"phylo")){
		stop("tree must be of class 'phylo'")
		}
	#simulate a perfect character dataset (parsimony informative binary chars) for a given tree
	charMat<-matrix(0,Ntip(tree),nchar)
	rownames(charMat)<-tree$tip.label
	desc<-sapply(prop.part(tree),function(x) tree$tip.label[x])
	desc<-desc[sapply(desc,length)!=Ntip(tree)]	        #get rid of root node that contains all taxa (not pars informative!
	nnode<-length(desc)
	if(nchar>nnode){
		#repeat desc if nchar multiple of nnode
		if((nchar %/% nnode) >1){
			for(i in 1:((nchar %/% nnode)-1) ){
				desc[(length(desc)+1):(length(desc)+nnode)]<-desc[1:nnode]
				}
			}
		if((nchar%%length(desc)) != 0){
			nrand<-nchar-length(desc)
			#sample without replacement
			desc[(length(desc)+1):nchar]<-desc[sample(1:nnode,nrand,replace=FALSE)]
			message(paste("Randomly sampling nodes for",nrand,"extra characters"))
			}
	}else{
		if(nnode>nchar){
			message("nchar less than the number of nodes, all returned characters placed randomly")
			#sample without replacement
			desc[1:nchar]<-desc[sample(1:nnode,nchar,replace=FALSE)]
			}
		}
	for(i in 1:nchar){
		charMat[desc[[i]],i]<-1
		}		
	return(charMat)
	}