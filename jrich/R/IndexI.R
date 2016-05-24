#'
#' @title I index value for a single topology.
#'
#' @description This function assigns the same weight to sister clades
#' (see Van-Wright et al., 1981). The input tree is reordered in post order.
#' 
#' @param tree is a single tree with n terminals, an ape phylo object.
#'
#' @return Returns a vector with weights.
#' 
#' @examples
#'  library(jrich)
#'  data(tree)
#'  plot(tree)
#'  indexi               <- IndexI(tree)
#'  newTree              <- tree
#'  newTree$tip.label    <- indexi
#'  plot(newTree)
#'  
#'
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'





IndexI <- function(tree=tree) {
  
  tree <- reorder.phylo(tree,order="postorder")
    
  matriz = matrix(0,nrow=1,ncol=(length(tree$tip.label)*2-1))
   

## masked functions


	Sisters <-
		function (tree=tree, node=node) {
  
  		# get the ancestor
  		ancestro<-tree$edge[(tree$edge[,2] == node)][1]
  
  		# get children of that ancestor
  		sisters<-tree$edge[c(which(tree$edge[,1] == ancestro)),2]
  
  	return(sisters)
	}



	Sisters.tip <-
		function (tree=tree, node=node) {
  		##
  		## recibe un arbol y un nodo y reporta si todos los hijos del ANCESTRO de ese nodo
  		## son tips
  		##  
  
  	if (all(Sisters(tree,node) <= length(tree$tip.label))) {
    	return(TRUE)
  		}else{
    return(FALSE)
  		}  
  
	}



	#

	Weight.sister.tips <-
	function (tree=tree, matriz=matriz) {
  
  		## llena de 1 la matriz de pesos iniciales 
  		## para terminales hermanas
  
  	for (node in 1:length(tree$tip.label)){
    	if ((node <= length(tree$tip.label)) & Sisters.tip(tree,node)){
      	matriz[1,node]  <- 1
   
      		#! ancestro           <-    tree$edge[(tree$edge[,2] == node)][1]
      		#! matriz[1,ancestro] <-    matriz[1,ancestro]+1
      		#! print(c(node,ancestro))
      
    	}
  	}
  	return(matriz)
	}


	#

	Weight.other.nodes <-
	function (tree=tree,matriz=matriz) {
  
  		## recorrer los tips para asignar pesos que no han sido asignados
  
  	for (i in ((length(tree$tip)+tree$Nnode):(length(tree$tip)+1))){
    
    	if((any(Children(tree=tree,node=i) <= length (tree$tip.label))) &
       	(any(Children(tree=tree,node=i) >  length (tree$tip.label)))){
      
      		matriz[1,(Children(tree=tree,node=i))] <- max(matriz[1,Children(tree=tree,node=i)]) 
      		matriz[1,i] <- sum(matriz[1,Children(tree=tree,node=i)])  
      	}else{
      		matriz[1,i] <- sum(matriz[1,Children(tree=tree,node=i)])  
    	} 
  	}
  
  	return(matriz[1,1:length(tree$tip)])
  
	}


  
#! pares de tips como 1
  matriz <- Weight.sister.tips(tree,matriz)
  
  #! completado
  matriz <- Weight.other.nodes(tree,matriz)
    #! revisar si es mejor escribir un data frame 
    #! con los nombres de las especies y sus pesos
  
  return(matriz)
}
