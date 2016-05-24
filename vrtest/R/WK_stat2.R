WK_stat2 <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    mq <- matrix(NA, nrow=length(kvec), ncol=1)
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    mq[i,] <- WK_stat1(y,k)
    }

return(max(abs(mq)))
}
