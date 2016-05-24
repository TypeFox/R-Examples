fuller1 <-
function(x,p)
{
x <- as.matrix(x)
n <- nrow(x)
y <- x[(p+1):n,1]
xmat <- as.matrix(x[p:(n-1),1])
if(p > 1) {
z <- x[2:n,1]-x[(1:n-1),1]
index1 <- (p-1):(n-2) 
for( i in 1:(p-1))
{xmat <- cbind(xmat,z[index1])
index1 <- index1-1}
}
b <- solve( t(xmat) %*% xmat) %*% t(xmat) %*% y
e <- y - xmat %*% b
s2 <- sum(e^2)/(n-p)
cov <- solve(t(xmat) %*% xmat)*s2
se <- sqrt(cov[1,1])
tau <- (b[1,1]-1)/se

return(list(coef=b,se=se,tau=tau))
}
