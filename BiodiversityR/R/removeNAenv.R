`removeNAenv` <-
function(x,variable) {
    subs <- is.na(x[,variable])
    subs <- (subs==F)
    x <- x[subs,,drop=F]
    for (i in 1:ncol(x)) {
        if (is.factor(x[,i])) {x[,i] <- factor(x[,i][drop=T])}
    }
    return(x)
}

