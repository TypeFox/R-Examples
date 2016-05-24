testMutationalPatternAll <- function (mtx, PRINT=FALSE) {
    mtx <- mtx+1
    ngenes <- ncol (mtx)
    gene1.ind <- unlist (mapply (seq.int, 1, 1:(ngenes-1)))
    gene2.ind <- unlist (mapply (rep, 2:ngenes, 1:(ngenes-1)))
    mtx.ind.tmp <- cbind (gene1.ind, gene2.ind)
    pvalues <- t(apply (mtx.ind.tmp, 1, testMutationalPatternPair, mtx=mtx, PRINT=PRINT))
    colnames (pvalues) <- c("pValueLL", "pValueGL", "pValueLG", "pValueGG", "pValueME")
    pvalues <- data.frame (gene1=colnames(mtx)[mtx.ind.tmp[,1]], gene2=colnames(mtx)[mtx.ind.tmp[,2]], pvalues)
    return (pvalues)
}

