Wright_stat <-
function(y,k) 
{
    y <- as.matrix(y)
    n <- nrow(y)
    ranking <- as.matrix(rank(y))
    r1 <- (ranking - 0.5*(n+1) )/sqrt((n-1)*(n+1)/12)
    r2 <- qnorm(ranking/(n+1))
    s <- sign(y)
    s[ s == 0] <- -1

    R1 <- stat(r1,k) 
    R2 <- stat(r2,k)
    S1 <- stat(s,k) 
    
return(list(WR1=R1,WR2=R2,WS1=S1))
}
