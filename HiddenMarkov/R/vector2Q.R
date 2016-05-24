vector2Q <- function(p){
    #   map vector to mmpp Q rate matrix
    m <- 0.5 + sqrt(1+4*length(p))/2
    p <- exp(p)
    Q <- matrix(0, ncol=m, nrow=m)
    k <- 1
    for (j in 1:m){
        for (i in 1:m){
            if (i!=j){
                Q[i,j] <- p[k]
                k <- k+1
            }
        }
    }
    diag(Q) <- -(Q %*% rep(1, m))
    return(Q)
}
