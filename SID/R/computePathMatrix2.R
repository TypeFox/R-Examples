computePathMatrix2 <- function(G,condSet,PathMatrix1, spars=FALSE)
# Copyright (c) 2013 - 2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    # The only difference to the function computePathMatrix is that this function changes
    # the graph by removing all edges that leave condSet.
    # If condSet is empty, it just returns PathMatrix1.
    
    p <- dim(G)[2]
    if(length(condSet) > 0)
    {
        G[condSet,]<-matrix(0,length(condSet),p)
        
        if(spars)
        {
            G <- Matrix(G)
            PathMatrix2 <- Diagonal(p) + G
        } else
        {
            PathMatrix2 <- diag(1,p) + G
            #PathMatrix2 <- as(diag(1,p) + G, "sparseMatrix")
            #   sparseMatrix does not seem to lead to a big improvement in speed (in fact, it seems slightly slower)
        }
        
        
        k <- ceiling(log(p)/log(2))
        for(i in 1:k)
        {
            PathMatrix2 <- PathMatrix2 %*% PathMatrix2            
        }
        PathMatrix2 <- PathMatrix2 > 0        
    } else
    {
        PathMatrix2 <- PathMatrix1
    }
    return(PathMatrix2)
}
