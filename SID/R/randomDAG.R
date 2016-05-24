randomDAG <-
function(p,probConnect,causalOrder = sample(p,p,replace=FALSE))
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
#    
# All rights reserved.  See the file COPYING for license terms. 
#
# simulates a directed acyclic graph (DAG) and returns its adjacency matrix
# 
# INPUT:
#   p           number of nodes 
#   probConnect the probability that an edge i -> j is added to the DAG
#   causalOrder starting with sink node (also called topological order)
#   
# OUTPUT:
#   DAG         Adjacency matrix of a directed acyclic graph (DAG)    
{
    #DAG <- matrix(0, p, p)
    DAG <- Matrix(0, p, p, sparse=TRUE)
    for(i in 1:(p-2))
    {
        node <- causalOrder[i]
        possibleParents <- causalOrder[(i+1):p]
        numberParents <- rbinom(n=1,size=(p-i),prob=probConnect)
        Parents <- sample(x = possibleParents, size = numberParents, replace = FALSE)
        DAG[Parents,node] <- rep(1,numberParents)
    }
    # Sample does not work properly when choosing from sets with one element. We thus consider the last case separately.  
    node <- causalOrder[p-1]
    ParentYesNo <- rbinom(n=1,size=1,prob=probConnect)
    DAG[causalOrder[p],node] <- ParentYesNo
    
    return(DAG)
}
