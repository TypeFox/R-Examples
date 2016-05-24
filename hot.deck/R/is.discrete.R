is.discrete <-
function(data, cutoff = 10){
    discrete <- vector()
    if(is.vector(data)){
    data <- na.omit(data)
    discrete <- length(unique(data)) <= cutoff | is.factor(data)
    }
    if(is.data.frame(data) | is.matrix(data)){
        for(j in 1:ncol(data)){
            tmp <- na.omit(data[,j])
            discrete[[j]] <- length(unique(tmp)) <= cutoff | is.factor(tmp)
        } 
    }
    return(discrete)
}
