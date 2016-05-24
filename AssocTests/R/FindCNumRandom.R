## Bootstrap to determine the number of groups
FindCNumRandom <- function(x, n, kG, n.monterCarlo)
{  
    # x             --> data matrix
    # n             --> nrow(x)
    # kG            --> number of total clusters
    # n.monterCarlo --> simulation times

    # centralize x
    x <- scale(x, center=TRUE, scale=FALSE)

    W <- CalculateWAll(x, n, kG)

    WStar <- matrix(data=0, nrow=n.monterCarlo, ncol=kG) 
    
    bound <- rbind(n, apply(x,2,range))

    for (i in 1:n.monterCarlo)
    {
        y <- apply(bound, 2, UniformSample)
        WStar[i,] <- CalculateWAll(y, n, kG)
    }

    CalculateGapK(W, WStar, kG, n.monterCarlo)
}
