Wald <-
function(y,kvec)
{
    y <- as.matrix(y)
    n <- nrow(y)
    mvr <- matrix(NA, nrow=length(kvec), ncol=1)
    for (i in 1:length(kvec))
    {
    k <- kvec[i]
    VR <- Wald1(y,k)
    mvr[i,] <- cbind(VR)
    }
    mat <-covmat(kvec)
    w <- n* t(mvr) %*% solve(mat) %*% mvr
    alpha <- c(0.1,0.05,0.01)
    cr <- qchisq(1-alpha,length(kvec))
    return(list(Holding.Period=kvec,Wald.stat=as.numeric(w),Critical.Values_10_5_1_percent=cr))
}
