thin<-function(x, k){
    ## returns every kth element of a vector, matrix, or data.frame

    if(is.vector(x)){
        n <- length(x)
        idx <- which((1:n)%%k==0)
        return(x[idx])
    }else if(is.matrix(x) | is.data.frame(x)){
        nRow <- nrow(x)
        idx <- which((1:nRow)%%k==0)
        return(x[idx,])
    }else{
        stop("x must be a vector, a matrix or a data.frame")
    }
}
