is.upset <-
function(z, Q=1) {
    if(is.numeric(Q)) {
        v <- 1:nrow(z)
        Q <- v %in% Q
    }
    if(is.character(Q))
        Q <- rownames(z) %in% Q
    all(Q == upset(z, Q))
}
