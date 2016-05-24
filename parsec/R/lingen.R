lingen <-
function(z) {
    res <- order(-levels.incidence(z))
    names(res) <- rownames(z)[res][res]
    return(res)
}
