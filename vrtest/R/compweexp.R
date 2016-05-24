compweexp <-
function(inf) 
{
n <- length(inf)
weiexp <- matrix(1,nrow=n,ncol=n)
for(i in 1:n)
{
    for(j in (i+1):n)
    {
    if(j > n) break
    aux1 <- (inf[i]-inf[j]) %*% t(inf[i]-inf[j])
    weiexp[i,j] <- exp(-0.5*aux1)
    weiexp[j,i] <- weiexp[i,j]
    }
}
return(weiexp)
}
