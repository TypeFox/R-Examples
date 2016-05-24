`allExpressions` <-
function(noflevels, raw=FALSE, arrange=FALSE) {
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    aEmat <- createMatrix(noflevels + 1)
    if (arrange) {
        aEmat <- QCA::sortMatrix(aEmat)
        sum.zeros <- apply(aEmat, 1, function(idx) sum(idx == 0))
        aEmat <- aEmat[order(sum.zeros, decreasing=TRUE), ]
    }
    return(structure(list(aE=aEmat - 1, raw=raw), class = "aE"))
}

