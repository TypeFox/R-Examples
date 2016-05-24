## __________________________________________________________
##
## BuildNetwork
##
## Given a list of scored edges and a REQUIRED vector of labels
## (e.g, 1:p), this function builds an object "Network". 
##
## This is useful for exporting a network described by a list of
## edges to a network described by an adjacency matrix.
##
##
## INPUT
##         - Edges : a n x 3 matrix (each line contains a couple of
##         vertices plus the associated score)
##         - Labels : a vector
##         - nonedges.val : optional. Value attributed to the not existing 
##         edges in the generated score matrix out\$Score, default=NA.
##
## OUTPUT
##         - Network : a list that contains the number of vertices, a
##         vector of labels of the vertices, a vector of the connected
##         vertices, the proportion of edges, the number of edges, the 
##         graph of the network (an adjacency matrix and a list of edges)
##         and a valued matrix for the edges.
## __________________________________________________________
##

BuildNetwork <- function(Edges,Labels,nonedges.val=NA){

  ## initialize
  Net = list()
  Net$Vertices = list()
  Net$Edges = list()
  
  ## vertices (genes)
  Net$Vertices$Num <- length(Labels)
  Net$Vertices$Labels <- Labels
  
  ## proportion and number of edges
  Net$Edges$Prop <- length(Edges[,1])/(length(Labels)^2)
  Net$Edges$Num <- dim(Edges)[1]
  
  ## matrix of edges (plus the associated weight)
  Net$Edges$List <- Edges

  ## indicates the "connected" genes of the network
  Net$Vertices$Connected <- which(Labels %in% Net$Edges$List[,1:2])

  ## Build the adjacency Matrix
  Net$AdjMatrix <- matrix(0,Net$Vertices$Num,Net$Vertices$Num)
  for (k in 1:Net$Edges$Num) {
    Net$AdjMatrix[which(Labels == Net$Edges$List[k,2]),which(Labels == Net$Edges$List[k,1])] <- 1
  }

  ## Build the matrix of  edge values Matrix
  Net$Score <- matrix(nonedges.val,Net$Vertices$Num,Net$Vertices$Num)
  for (k in 1:Net$Edges$Num) {
    Net$Score[which(Labels == Net$Edges$List[k,2]),which(Labels == Net$Edges$List[k,1])] <- Edges[k,3]
  }
  
  return(Net)
}
