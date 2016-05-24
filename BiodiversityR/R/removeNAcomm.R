`removeNAcomm` <-
function(x,y,variable) {
    subs <- is.na(y[,variable])
    subs <- (subs==F)
    x <- x[subs,,drop=F]
    return(x)
}


