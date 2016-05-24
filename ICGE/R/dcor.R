dcor<- function(x){
########## Correlation distance between pais of objects ##################
# Input:
# x: data matrix
#
# Output:
# d: distance matrix
###########################################################################
    n<- dim(x)[1]
    p <- dim(x)[2]

    d<- matrix(0, n,n)
    for (i in 1:(n-1)){
         for (j in (i+1):n){
            r<- cor(x[i,],x[j,])
            d[i,j] <- sqrt(1-r)
            d[j,i] <- d[i,j]
         }
     }

    d <- as.dist(d)
    return(d)
}
