transitiveClosure.cover <-
function(m) {
    res <- transitiveClosure.default(m)
    class(res) <- "cover"
    return(res)
}
