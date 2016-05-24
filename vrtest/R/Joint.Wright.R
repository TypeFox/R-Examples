Joint.Wright <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    W_mat <- matrix(NA, nrow=length(kvec), ncol=3)
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    W <- Wright_stat(y,k)
    W_mat[i,] <- cbind(W$WR1,W$WR2,W$WS1)
    }

    jr1 <- max(abs(W_mat[,1]))
    jr2 <- max(abs(W_mat[,2]))
    js1 <- max(abs(W_mat[,3]))
    
return(list(Holding.Period=kvec,JR1=jr1,JR2=jr2,JS1=js1))
}
