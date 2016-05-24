## Author: Jim Young

## Version of 'makesimprob' that samples from Dirichlet posterior rather than

# using the expected value - see 'Bayesian data analysis' Gelman,
# Carlin, Stern and Rubin 1995 p482.

fullsimprob <- function (nw) {
    for (nid in 1:nw$n) {
        node <- nw$nodes[[nid]]
        if (node$type == "continuous") stop("fullsimprob only works for discrete
nodes")
        parents <- node$parents
        dparents <- c()
        if (nw$nd > 0) dparents <- sort(intersect(parents, nw$discrete))
        if (length(dparents) > 0) {
            Dim <- c()
            dnames <- list(node$levelnames)
            for (i in dparents) {
                Dim <- c(Dim, nw$nodes[[i]]$levels)
                dnames <- c(dnames, list(nw$nodes[[i]]$levelnames))
            }
        }
        if (identical(length(dparents),as.integer(0))) {
            dnames <- list(node$levelnames)
            Dim <- c()
        }
# Additional code to extract conditional posterior frequencies and re-sort.
# Set up an empty array to hold dimensions of conditional posterior.
        CDim <- c()
# Re-order node and its discrete parents into network order.
        netorder <- sort(union(nid, dparents))
# Find dimensions of these nodes in network order.
        for (i in netorder) {
             CDim <- c(CDim, nw$nodes[[i]]$levels)
        }
# Pull conditonal posterior counts out of the network for this node.
        condP <- array(unlist(node$condposterior),dim=CDim)
# Sampling from a Dirichlet distribution with parameters alpha - see Gelman
# Carlin, Stern and Rubin p 482: draw x's from independent gamma distributions
# with shape parameters alpha and common scale, then thetas equal each x

# divided by the sum of all x's.
        condP <- array(rgamma(n=length(unlist(condP)),shape=condP,scale=1),
                       dim=CDim)
        condP <- condP/sum(unlist(condP))
# End of Dirichlet sampling code. The last line appears to be unnecessary.
# The next line is critical ? reorder into dimensions expected by ?rnetwork?.
# If there are parents, turn array around using order of node then parents.
        if (length(dparents) > 0) condP <- aperm(condP, rank(c(nid,dparents)))
        Dim <- c(node$levels, Dim)
#  Next line, instead of ?simtab <- array(1/prod(Dim), dim = Dim)?
#  use the conditional posterior probabilities?
        simtab <- condP
        dimnames(simtab) <- dnames
        if (length(node$parents) > 0) 
            simtab <- prop.table(simtab, 2:(length(node$parents) + 1))
            else simtab <- prop.table(simtab)
        nw$nodes[[nid]]$simprob <- simtab
     }
    nw
}


