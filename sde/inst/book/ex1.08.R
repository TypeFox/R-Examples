# ex1.08.R
set.seed(123)
phi <- function(i,t,T){
   (2*sqrt(2*T))/((2*i+1)*pi) * sin(((2*i+1)*pi*t)/(2*T))
}
T <- 1
N <- 100
t <- seq(0,T,length=N+1)
W <- numeric(N+1)
n <- 10
Z <- rnorm(n)
for(i in 2:(N+1)) 
   W[i] <- sum(Z*sapply(1:n, function(x) phi(x,t[i],T)))
plot(t,W,type="l",ylim=c(-1,1))
n <- 50
Z <- rnorm(n)
for(i in 2:(N+1)) 
   W[i] <- sum(Z*sapply(1:n, function(x) phi(x,t[i],T)))
lines(t,W,lty=2)
n <- 100
Z <- rnorm(n)
for(i in 2:(N+1)) 
   W[i] <- sum(Z*sapply(1:n, function(x) phi(x,t[i],T)))
lines(t,W,lty=3)
