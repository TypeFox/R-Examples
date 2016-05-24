#'
#' @title W index value for a single topology.
#'
#' @description This function assigns the weight according to the ramification patterns 
#' (see Van-Wright et al., 1981). 
#' The input tree is reordered in post order.
#' Returns a vector with weights.
#' 
#' @param tree is a single tree with n terminals, an ape phylo object.
#'
#' @examples
#'   library(jrich)
#'   data(tree)
#'   plot(tree)
#'   indexw             <- IndexW(tree)
#'   newTree            <- tree
#'   newTree$tip.label  <- indexw
#'   plot(newTree)
#'   
#'
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'





IndexW <- function (tree=tree) {

	## masked functions


	Ancestor <- function (tree=tree,node=node) {
  
  		ancestro <- tree$edge[(tree$edge[,2] == node)][1]
  
  	return(ancestro)
  
	}


  
  node.Matrix  <-  matrix(0,nrow=1,ncol=(length(tree$tip.label)*2-1))
  
  tree <- reorder.phylo(tree,order="postorder")
  
  raiz <-length(tree$tip.label)+1
  
  node.Matrix[1,raiz] <- 1
  
  hijos        <- Children(tree=tree,node=raiz)
  
  node.Matrix[1,hijos] <- node.Matrix[1,raiz]
  
  for (i in (((length(tree$tip.label)+2)):(length(tree$tip.label)+tree$Nnode))){

    node.Matrix[1,i]     <-  node.Matrix[1,Ancestor(tree=tree,node=i)]+1 
    hijos           <-  Children(tree=tree,node=i)
    node.Matrix[1,hijos] <-  node.Matrix[1,i]
        
  }
  
  node.Matrix <- node.Matrix[1,1:length(tree$tip.label)]
  node.Matrix <- sum(node.Matrix)/node.Matrix
  node.Matrix <- node.Matrix/min(node.Matrix)

  return(node.Matrix)
}
