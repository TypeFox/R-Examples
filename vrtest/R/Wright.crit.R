Wright.crit <-
function(n,k,nit)
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

    R1 <- stat(r1,k) 
    R2 <- stat(r2,k)
    S1 <- stat(s,k) 
    mat[i,] <- c(R1,R2,S1)
    }

alpha <- c(0.01,0.05,0.1)
R1crit <- quantile(mat[,1],c(0.5*alpha, rev(1-0.5*alpha)) )
R2crit <- quantile(mat[,2],c(0.5*alpha, rev(1-0.5*alpha)) )
S1crit <- quantile(mat[,3],c(0.5*alpha, rev(1-0.5*alpha)) )
return(list(Holding.Period=k,R1.crit=R1crit,R2.crit=R2crit,S1.crit=S1crit))    
}
