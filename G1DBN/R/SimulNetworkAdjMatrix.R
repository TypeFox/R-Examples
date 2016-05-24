
## __________________________________________________________
##
## FUNCTION SimulNetworkAdjMatrix
##
## This function builds a object "Network" by simulating a
## matrix of valued adjacencies from a number of vertices, a proportion
## of edges and the range of the uniform distribution that is used to build
## the adjacency matrix. An optional vector of labels may be given.
##
## INPUT
##         - Num : number of genes
##         - EdgesProp : edges proportion in the network
##         - Range : vector with 4 elements specifying range values for the 
##         adjacency matrix generation (minimum negative value, maximum negative
##         value, minimum positive value, maximum positive value)
##         - Labels : an optional vector of labels for the edges
##
## OUTPUT
##         - Network : a list that contains
##            the number of vertices, a vector of labels
##            of the vertices, a vector of the regulated
##            vertices, the proportion of edges, the number 
##            of edges, an adjacency matrix (binary) 
##            and a valued adjacency matrix
## __________________________________________________________
##

SimulNetworkAdjMatrix <- function(Num,EdgesProp,Range,Labels=1:Num){
  
  ## initialize
  Net = list()
  Net$Vertices = list()
  Net$Edges = list()
  
  ## vertices (genes)
  Net$Vertices$Num <- Num
  Net$Vertices$Labels <- Labels
  
  ## proportion and number of of edges
  Net$Edges$Prop <- EdgesProp
  Net$Edges$Num <- floor(Num^2*EdgesProp)
  
  ##  generating matrix A
  m <- expand.grid(1:Num,1:Num)
  
  A <- matrix(0,Num,Num)
  for (k in sample(1:Num^2,Net$Edges$Num)) {
    A[m[k,1],m[k,2]] <- sample(c(runif(1,Range[1],Range[2]),
                                 runif(1,Range[3],Range[4])),1)
  }
  Net$A=A

  ##ADJACENCY MATRIX
  Net$AdjMatrix = (A!=0)*1



  ## indicates the regulated genes of the network
  Net$Vertices$Regulated <- which(A %*% matrix(1,Num,1)
               + t(A) %*% matrix(1,Num,1) != 0)

  return(Net)
}
