`same.sites` <-
function(x,y) {    
    n <- nrow(y)
    p <- ncol(x)
    result <- array(0,dim=c(n,p))
    for (i in 1:n) {
        index <- rownames(y)[i] == rownames(x)
        if (any(index)==T) {
            result[i,] <- as.matrix(x[index,,drop=F])
        }
    }  
    result <- data.frame(result)
    rownames(result) <- rownames(y)
    colnames(result) <- colnames(x)
    totals1 <- rowSums(result)
    index1 <- totals1 == 0
    zerosites1 <- rownames(result)[index1]
    totals2 <- rowSums(x)
    index2 <- totals2 == 0
    zerosites2 <- rownames(x)[index2]
    if (length(zerosites1)>0 || length(zerosites2) >0) {
        if (any(zerosites1==zerosites2)==F) {
            cat("Warning: some sites without species are different in original and resulting data\n")
            cat("Original sites without species: ", zerosites2, "\n")
            cat("Resulting sites without species: ", zerosites1, "\n")        
        }
    }
    return(result)
}

