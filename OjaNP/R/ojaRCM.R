`ojaRCM` <-
function(X, p = NULL, silent = FALSE, na.action = na.fail){
    OjaRanks <- ojaRank(X, x = NULL, p = p, silent = silent, na.action = na.fail)
    RCM <- crossprod(OjaRanks)/nrow(OjaRanks)
    dimnames(RCM) <- NULL
    rownames(RCM) <- colnames(X)
    colnames(RCM) <- colnames(X)
    return(RCM)
}
