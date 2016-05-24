JWright.crit <-
function(n,kvec,nit)
{
    set.seed(12345)
    mat <- matrix(NA,nrow=nit,ncol=3)
    for (i in 1:nit)
    {
        
        ranking <- as.matrix(sample(1:n,n,replace=FALSE))
        r1 <- (ranking - 0.5*(n+1) )/sqrt((n-1)*(n+1)/12)
        r2 <- qnorm(ranking/(n+1))
        y <- as.matrix(rnorm(n))
        s <- sign(y)
        s[ s == 0] <- -1
    
        statmat <- matrix(NA, nrow=length(kvec), ncol=3)
    
        for (j in 1:length(kvec))
        {
        k <- kvec[j]
        statmat[j,] <- cbind(stat(r1,k),stat(r2,k),stat(s,k))
        }
    
    R1 <- max(abs(statmat[,1]))
    R2 <- max(abs(statmat[,2]))
    S1 <- max(abs(statmat[,3]))
    mat[i,] <- c(R1,R2,S1)
    }

alpha <- c(0.90,0.95,0.99)
R1crit <- quantile(mat[,1],alpha )
R2crit <- quantile(mat[,2],alpha )
S1crit <- quantile(mat[,3],alpha )
return(list(Holding.Period=kvec,JR1.crit=R1crit,JR2.crit=R2crit,JS1.crit=S1crit))    
}
