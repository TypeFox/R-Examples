#' Compute small-worldness relative to random graphs.
#' 
#' Function generates random graphs for the same number of nodes and edges, and
#' estimates small-worldness.
#'
#' @param g A graph object.
#' @return a design matrix.
#' @export
#'

computeSmallWorldness <- function (g, n, m, nrep=1000) {
    g.rand.clust <- rep (0, nrep)
    g.rand.l <- rep (0, nrep)
    for (i in 1:nrep) {
        g.rand.tmp <- erdos.renyi.game (n, m, type="gnm")
        g.rand.clust[i] <- transitivity (g.rand.tmp)
        g.rand.l[i] <- average.path.length (g.rand.tmp)
    }
    
    return (transitivity (g) / mean (g.rand.clust) * mean (g.rand.l) / average.path.length (g))
    
}

