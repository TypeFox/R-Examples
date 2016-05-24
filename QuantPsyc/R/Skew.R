"Skew" <-
function(x)  
{
n <- length (x[!(is.na(x))])
sd <- sqrt(var(x,na.rm=TRUE))
m <- mean(x,na.rm=TRUE)
sk <- (n/((n-1)*(n-2)))*sum(((x-m)/sd)^3, na.rm=TRUE)
se.sk <- sqrt(6/n)
t.sk <- sk/se.sk
p.sk <- 1-pnorm(abs(t.sk))#Note - evaluated against z instead of t
matSk <- cbind(sk,se.sk,t.sk,p.sk)
return(matSk)
}

