vector2Pi <- function(p){
    #   map vector to dthmm transition prob matrix
    m <- 0.5 + sqrt(1+4*length(p))/2
    invlogit <- function(eta)
        exp(eta)/(1+exp(eta))
    #   last column of Pi constrained by others
    Pi <- matrix(0, ncol=m, nrow=m)
    k <- 1
    for (i in 1:m){
        for (j in 1:(m-1)){
            Pi[i,j] <- (1-sum(Pi[i,]))*invlogit(p[k])
            k <- k+1
        }
        Pi[i,m] <- 1 - sum(Pi[i,])
    }
    return(Pi)
}
