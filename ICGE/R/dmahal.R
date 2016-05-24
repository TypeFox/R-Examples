dmahal<- function(datos, S){
########## Mahalanobis distance between pais of objects ##################
# Input:
# datos: data matrix
# S: var-cov matrix
# Output:
# d: distance matrix
###########################################################################
    n<- dim(datos)[1]
    p <- dim(datos)[2]
    s1 <- dim(S)[1]
    s2 <- dim(S)[2]
    if (s1 != s2){stop("S must be squared matrix")}
    if(p != s1){stop("data and covariance matrix must have coherent dimesions")}
    d<- matrix(0, n,n)
    I <- diag(rep(1,p)) #identity matrix 
    Sinv <- solve(S, I)
    for (i in 1:n){
         for (j in 1:i){
            a<- datos[i,]-datos[j,]
            a<-matrix(a, ncol=1)
            d[i,j] <- sqrt(t(a)%*%Sinv%*%a)
            d[j,i] <- d[i,j]
         }
     }

d <- as.dist(d)
return(d)
}
