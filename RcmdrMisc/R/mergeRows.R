# simple row-wise merge of data frames

# last modified 2014-08-04 by J. Fox

mergeRows <- function(X, Y, common.only=FALSE, ...){
    UseMethod("mergeRows")
}

mergeRows.data.frame <- function(X, Y, common.only=FALSE, ...){
    cols1 <- names(X)
    cols2 <- names(Y)
    if (common.only){
        common <- intersect(cols1, cols2)
        rbind(X[, common], Y[, common])
    }
    else {
        all <- union(cols1, cols2)
        miss1 <- setdiff(all, cols1)
        miss2 <- setdiff(all, cols2)
        X[, miss1] <- NA
        Y[, miss2] <- NA
        rbind(X, Y)
    }
}