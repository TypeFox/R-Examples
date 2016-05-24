hipamAnthropom <- function(data,asw.tol=0,maxsplit=5,local.const=NULL,orness=0.7,type,
                           ah=c(23,28,20,25,25),verbose,...){
 #Initialize the tree:
 tree <- initialize.tree(data, maxsplit, orness, type, ah, verbose, ...)
 #Local hipam: 
 tree <- hipam.local(tree, data, asw.tol, maxsplit, local.const, orness, type, ah, verbose, ...)
 tree$cases <- tree$medoids
 tree$medoids <- NULL
 tree
}


