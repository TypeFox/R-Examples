##
## Function for creating 'NNOC' objects 
nnoc <- function(G, h){
    return(list(conType = "NNOC", G = as.matrix(G), h = as.matrix(h), dims = nrow(G)))
}
