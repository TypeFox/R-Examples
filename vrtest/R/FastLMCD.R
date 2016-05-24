FastLMCD <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    mq <- matrix(NA, nrow=length(kvec), ncol=2)
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    LM <- FastLM_stat(y,k)
    mq[i,] <- cbind(LM$LM1,LM$LM2)
    }

    mv1 <- max(abs(mq[,1]))
    mv2 <- max(abs(mq[,2]))
    
return(list(M1=mq[,1],M2=mq[,2],CD1=mv1,CD2=mv2))
}
