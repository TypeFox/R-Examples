VAR1.sim <-
function(a0,a1,n,sigu)
{
nob <- 50+n
x <- matrix(0,nrow=nob,ncol=2)
tem <- chol(sigu)
for (i in 2:nob)
x[i,] <-  a0 + a1 %*% matrix(x[i-1,]) + crossprod(tem, matrix(rnorm(2),nrow=2))
return(x[51:nob,])
}
