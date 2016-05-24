testMutationalPatternAll.wrapper <- function (data, QVALUE=TRUE, PRINT=FALSE) {
    mtx <- as.matrix (t(data[,-1]))
    colnames (mtx) <- data[,1]

    pvalues <- testMutationalPatternAll (mtx, PRINT=PRINT)
    if (QVALUE) {
        qvalues <- matrix (qvalue (as.vector (as.matrix(pvalues[,3:7])))$qvalue, byrow=FALSE, ncol=5)
        colnames (qvalues) <- c("qValueLL", "qValueGL", "qValueLG", "qValueGG", "qValueME")
        qvalues <- data.frame (pvalues[,1:2], qvalues)
        return (list (pvalues = pvalues, qvalues = qvalues))
    }
    else {
        return (list (pvalues=pvalues))
    }
}

