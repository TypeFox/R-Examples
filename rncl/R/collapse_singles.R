## This is the collapse.singles function from ape 3.3
## written by Emmanuel Paradis
## I included here a copy of it as the current testing version of ape 3.3.0.5
## has a bug in it. This version seems to work for rncl at this stage.
collapse_singles <- function(tree) {
    elen <- tree$edge.length
    xmat <- tree$edge
    node.lab <- tree$node.label
    nnode <- tree$Nnode
    ntip <- length(tree$tip.label)
    singles <- NA
    while (length(singles) > 0) {
        tx <- tabulate(xmat[, 1])
        singles <- which(tx == 1)
        if (length(singles) > 0) {
            i <- singles[1]
            prev.node <- which(xmat[, 2] == i)
            next.node <- which(xmat[, 1] == i)
            xmat[prev.node, 2] <- xmat[next.node, 2]
            xmat <- xmat[xmat[, 1] != i, ]
            xmat[xmat > i] <- xmat[xmat > i] - 1L
            elen[prev.node] <- elen[prev.node] + elen[next.node]
            if (!is.null(node.lab))
                node.lab <- node.lab[-c(i - ntip)]
            nnode <- nnode - 1L
            elen <- elen[-next.node]
        }
    }
    tree$edge <- xmat
    tree$edge.length <- elen
    tree$node.label <- node.lab
    tree$Nnode <- nnode
    tree
}
