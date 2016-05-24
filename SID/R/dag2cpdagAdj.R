dag2cpdagAdj <-
function(Adj)
# Copyright (c) 2010 - 2012  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
{
    if(sum(Adj) == 0)
    {
        return(Adj)
    }
    cO <- computeCausOrder(Adj)
    d <- as(Adj[cO,cO], "graphNEL")
    cpd <- pcalg::dag2cpdag(d)
    res <- matrix(NA,dim(Adj)[1],dim(Adj)[1])
    res[cO,cO] <- as(cpd, "matrix")
    result <- res
    ################
    #THE CODE ABOVE USES THE CAUSAL ORDER BECAUSE OF A VERY WEIRD BEHAVIOUR IN PCALG!!!
    ################
    ##Adj <- cbind(c(0,0,0,0),c(1,0,0,0),c(1,1,0,0),c(1,1,1,0))
    ##as(pcalg::dag2cpdag(as(Adj, "graphNEL")), "matrix")
    ##plot(pcalg::dag2cpdag(as(Adj, "graphNEL")))
    ##as(pcalg::dag2cpdag(as(t(Adj), "graphNEL")), "matrix")
    ##plot(pcalg::dag2cpdag(as(t(Adj), "graphNEL")))
    
    
    
    return(result)
}
