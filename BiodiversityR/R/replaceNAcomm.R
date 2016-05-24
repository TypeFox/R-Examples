`replaceNAcomm` <-
function(x) {
    for (j in 1:ncol(x)) {
        x[is.na(x[,j]), j] <- 0 
    }
    return(x)
}


