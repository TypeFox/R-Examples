Lo.Mac <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    mq <- matrix(NA, nrow=length(kvec), ncol=2)
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    LM <- LM_stat(y,k)
    mq[i,] <- cbind(LM$LM1,LM$LM2)
    }
    VR <- mq
    rownames(VR) <- paste("k=",kvec,sep="")
    colnames(VR) <- c("M1","M2")
    return(list(Stats=VR))
    #return(list(Holding.Periods=kvec,M1.stat=mq[,1],M2.stat=mq[,2]))
}
