alphas0 <-
function(alpha,b,p,n)
{set.seed(12345)
b[1,1] <- alpha
a <- arlevel(b,p)

stat <- matrix(NA,nrow=500)
for(ii in 1:500)
{
nob <- n+50
u <- rnorm(nob)
y <- matrix(0,nrow=p)
    for( i in (p+1):nob)
    {tem <- sum( a[1:p,1] * y[(i-1):(i-p),1] ) + u[i]
    y <- rbind(y,tem)}
xs <- y[(nob-n+1):nob,1]
bs <- LSE(xs,p)$coef
stat[ii,1] <- bs[1]
}
return(median(stat))
}
