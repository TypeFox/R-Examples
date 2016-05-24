## calculate the distance between inner clusters, W
CalculateWAll <- function(x, n, kG)
{    
    W <- rep(1, kG)
     
    for (k in 2:kG)
    { 
        W[k] <- cluster::clara(x, k, samples=20, medoids.x=FALSE)$objective
    }

    W
}