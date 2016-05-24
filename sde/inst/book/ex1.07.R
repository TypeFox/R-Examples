# ex1.07.R
set.seed(123)
n <- 10 # far from the CLT
T <- 1
t <- seq(0,T,length=100)
S <- cumsum(2*(runif(n)>0.5)-1)	
W <- sapply(t, function(x) ifelse(n*x>0,S[n*x],0))
W <- as.numeric(W)/sqrt(n)
plot(t,W,type="l",ylim=c(-1,1))
n <- 100 # closer to the CLT
S <- cumsum(2*(runif(n)>0.5)-1)	
W <- sapply(t, function(x) ifelse(n*x>0,S[n*x],0))
W <- as.numeric(W)/sqrt(n)
lines(t,W,lty=2)
n <- 1000 # quite close to the limit
S <- cumsum(2*(runif(n)>0.5)-1)	
W <- sapply(t, function(x) ifelse(n*x>0,S[n*x],0))
W <- as.numeric(W)/sqrt(n)
lines(t,W,lty=3)
