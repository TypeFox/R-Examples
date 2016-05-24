##
## Function for creating 'NLFC' objects 
nlfc <- function(G, h){
    return(list(conType = "NLFC", G = as.matrix(G), h = as.matrix(h), dims = nrow(G)))
}
