moralize <- function (g)
{
    p <- vcount(g)
    g.moral <- as.undirected(g)
    A <- get.adjacency(g)
    diag(A) <- 0
    M <- A > 0
    Shared <- apply(M, 2, which)
    for (i in 1:p) {
        if (length(Shared[[i]]) > 1) 
            g.moral <- add.edges(g.moral, c(combn(Shared[[i]], 
                2)))
    }
    g.moral <- simplify(g.moral)
    return(g.moral)
}
