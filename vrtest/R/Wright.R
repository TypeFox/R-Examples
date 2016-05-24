Wright <-
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

    VR <- W_mat
    rownames(VR) <- paste("k=",kvec,sep="")
    colnames(VR) <- c("R1","R2","R3")
    return(list(Stats=VR))
}
