#' as.bipartite  --- convert a network object to a matrix for 
#' analysis with the bipartite package
#' INPUT = network model
#' OUTPUT = matrix representation
#' M. Lau July 2015
#' ---------------------------------------------------

as.bipartite <- function(x,y){
    y <- factor(y)
    unpack(x)$F[y == levels(y)[1],y == levels(y)[2]]
}

