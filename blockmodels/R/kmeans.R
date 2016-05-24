
blockmodelskmeans <- function(coordinates,k)
{

    n <- dim(coordinates)[1]
    dim <- dim(coordinates)[2]

    if(k==1)
    {
        return(c(matrix(1,n,1)))
    }

    centroids <- matrix(0,k,dim)

    dists <- as.matrix( dist( coordinates ))

    choosen <- which(dists==max(dists),arr.ind = TRUE)[1,c('row','col')]

    while(length(choosen)<k)
    {
        choosen<-c(choosen,which.max(apply(dists[choosen,],2,min)))
    }

    centroids <- coordinates[choosen,]

    #kmeans(coordinates,centroids,iter.max=n)$cluster;
    classif <- .Call("kmeans",
                     coordinates,
                     centroids,
                     PACKAGE="blockmodels")
    return(as.vector(1+classif))

}

