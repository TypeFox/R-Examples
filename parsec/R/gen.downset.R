gen.downset <-
function(z, Q=1) {
    if(is.numeric(Q)) {
        v <- 1:nrow(z)
        Q <- v %in% Q
    }
    if(is.character(Q))
        Q <- rownames(z) %in% Q
    if(!is.downset(z, Q))
        stop(paste(toString(rownames(z)[Q]), "is not a downset"))
    sub <- z[Q, Q]
    class(sub) <- class(z)
    Q[Q] <- maximal(sub)
    names(Q) <- rownames(z)
    return(Q)
}
