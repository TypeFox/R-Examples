`dist.zeroes` <-
function(comm,dist) {
    tots <- rowSums(comm)
    dist <- as.matrix(dist)
    l <- length(tots)
    for (i in 1:(l-1)) {
        if (tots[i]==0) {
            for (j in 2:l) {
                if ((tots[j]==0) && (is.finite(dist[i,j])==F)) {
                    dist[i,j] <- 0
                    dist[j,i] <- 0
                }
                if ((tots[j]>0) && (is.finite(dist[i,j])==F)) {
                    dist[i,j] <- 1
                    dist[j,i] <- 1
                }
             }       
        }
    }
    dist <- as.dist(dist)
    return(dist)
}

