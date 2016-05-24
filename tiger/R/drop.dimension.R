drop.dimension <- function(x){

    if(length(dim(x))==2){
        result <- vector(mode="numeric", length=dim(x)[1]*dim(x)[2])
        for(counter in 1:dim(x)[1]){
            start <- (counter-1)*dim(x)[2]+1
            end <- (counter)*dim(x)[2]
            result[start:end] <- x[counter,]
        }
        return(result)
    } 
    if(length(dim(x))==3){
        result <- array(dim=c(dim(x)[1]*dim(x)[2],dim(x)[3]))
        for(counter in 1:dim(x)[1]){
            start <- (counter-1)*dim(x)[2]+1
            end <- (counter)*dim(x)[2]
            result[start:end,] <- x[counter,,]
        }
        return(result)
    }
}
