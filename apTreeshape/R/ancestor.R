
ancestor <- function(tree, i){
which( (tree$merge[,1] == i) | (tree$merge[,2] == i) ) 
}

