#'
#' @title Children of a node.
#'
#' @description Get the children of a node in a tree.
#'
#' @param tree is a single tree with n terminals, an ape phylo object.
#' 
#' @param node representing the node in APE notation,is an integer.
#'
#' @return The children nodes of the internal node; in most cases, two integers.
#'
#'
#' @examples
#'  library(jrich)
#'  data(tree)
#'
#' Children(tree,7)
#'
#'@author Miranda-Esquivel Daniel R.
#'
#'






Children <-
function (tree=tree,node=node) {

  child <- tree$edge[c(which(tree$edge[,1] == node)),2]
  
  return(child)
}
