Chow.Denning <-
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

    mv1 <- max(abs(mq[,1]))
    mv2 <- max(abs(mq[,2]))
    
    alpha <- c(0.1,0.05,0.01)
    per <- 0.5*( 1-(1-alpha)^(1/length(kvec)))
    crit <- qnorm(1-per)

return(list(Holding.Periods=kvec,CD1=mv1,CD2=mv2,Critical.Values_10_5_1_percent=crit))
}
