transitiveClosure.default <-
function(m) {
    res <- m
    gk <- m
    n <- dim(m)[1]
    for(k in 2:n) {
        gk <- gk %*% m
        res <- res + gk
    }
    res <- res>0
    return(res)
}
