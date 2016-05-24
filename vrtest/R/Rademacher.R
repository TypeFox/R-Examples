Rademacher <-
function(n)
{
p <- 0.5
zmat <- rep(1,n);
u <- runif(n,0,1)
zmat[u > p] <- -1
return(zmat)
}
