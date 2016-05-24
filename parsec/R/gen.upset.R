gen.upset <-
function(z, Q=1) {
    if(is.numeric(Q)) {
        v <- 1:nrow(z)
        Q <- v %in% Q
    }
    if(is.character(Q))
        Q <- rownames(z) %in% Q
    if(!is.upset(z, Q)) stop(paste(toString(rownames(z)[Q]), "is not an upset"))
    sub <- z[Q, Q]
    class(sub) <- class(z)
    Q[Q] <- minimal(sub)
    names(Q) <- rownames(z)
    return(Q)
}
