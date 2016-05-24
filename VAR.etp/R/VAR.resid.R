VAR.resid <-
function(x,b,z,p)
{
n <- nrow(x)
k <- ncol(x)
y <- t(x[(p+1):n,])
e <- t(y - b%*%z)
tem <- colMeans(e)
es <- matrix(NA,nrow=nrow(e),ncol=ncol(e))
for( i in 1:k)
{
es[,i] <- matrix(e[,i] - mean(e[,i]))
}
return(es)
}
