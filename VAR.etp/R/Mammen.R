Mammen <-
function(n)
{
p <- (sqrt(5)+1)/(2*sqrt(5))
zmat <- rep(1,n)*(-(sqrt(5)-1)/2);
u <- runif(n,0,1)
zmat[u > p] <- (sqrt(5)+1)/2
return(zmat)
}
