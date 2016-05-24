ikcirt.data2Y <-
function(mxData, mxDelta) {
    
    mx.out2 <- mxData ; mx.out2[ is.na(mx.out2) ] <- 2
    
    df <- mx.out2 ; df
    
    Y <- mxDelta %*% df ; Y
    
    zeromask <- Y == 0
    Y[Y > 0.5] <- 1
    Y[Y < -0.5] <- 0
    Y[zeromask] <- NA
    
    mxData
    Y
    return(Y)
}
