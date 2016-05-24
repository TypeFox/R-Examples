dp.post.est <- 
function(x,y,alpha, lambda=mean(y))
{
n<-length(y)
c1 <- alpha/(alpha+n) 
c2 <- n/(alpha+n)
post.pmf <- rep(0, length(x))
for(i in 1:length(x)) {
post.pmf[i] <- c1* dpois(x[i],lambda) + c2 * (1/n) * sum(y==x[i])
}
post.pmf
}

