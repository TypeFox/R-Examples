getzeta.wprof <-
function(y) {
    y <- y$profiles
    m <- nrow(y)
    nam <- rownames(y)
    res <- matrix(NA, m, m, dimnames=list(nam, nam))
    for(a in 1:m) for(b in 1:m)
        res[a, b] <- all(y[a,]<=y[b,])
    class(res) <- "incidence"
    return(res)
}
