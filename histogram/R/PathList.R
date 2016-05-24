`PathList` <-
function(C,D){
    L <- nrow(C)
    for(i in D:1) {
        L <- c(C[L[1],i],L)
    }
    return(L)
}

