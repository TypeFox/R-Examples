Pi2vector <- function(Pi){
    #   map dthmm transition prob matrix to vector
    #   use logit like constraints
    m <- nrow(Pi)
    p <- rep(0, m*(m-1))
    #   last column of Pi constrained by others
    k <- 1
    for (i in 1:m){
        for (j in 1:(m-1)){
            p[k] <- log(Pi[i,j]/(1-sum(Pi[i,(1:j)])))
            k <- k+1
        }
    }
    return(p)
}
