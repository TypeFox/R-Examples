`removezerospecies` <-
function(x) {
    freq <- apply(x,2,sum)
    subs <- freq>0
    x <- x[,subs,drop=F]
    return(x)
}
