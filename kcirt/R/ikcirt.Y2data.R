ikcirt.Y2data <-
function(Y, mxDelta, ns) {
    
    Z <- 2*Y - 1
    Z[is.na(Z)] <- 2000
    
    df2 <- t(mxDelta) %*% Z
    df2
    
    
    mxData <- NULL
    ndxVec <- c(0, cumsum(ns))+1 ; ndxVec
    
    for(i in 1:(length(ns))) {
        
        dfa <- df2[ ndxVec[i]:(ndxVec[i+1]-1) ,  ]
        dfb <- dfa
        dfb[ dfb < -1000 | dfb > 1000 ] <- 0
        mxi <- apply( dfb, 2, function(x) { return(rank(x, na.last=TRUE)) } ) ; mxi
        mxi[ dfa < -1000 | dfa > 1000 ] <- NA ; mxi
        mxData <- rbind(mxData, mxi)
    }
    
    mxData
    return(mxData)
}
