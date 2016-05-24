# ex1.09.R
set.seed(123)
phi <- function(i,t,T){
   (2*sqrt(2*T))/((2*i+1)*pi) * sin(((2*i+1)*pi*t)/(2*T))
}
n <- 100
Z <- rnorm(n)
Delta <- seq(1e-7, 1e-2,length=30)
W <- sum(Z*sapply(1:n, function(x) phi(x,0.5,T)))
for(i in Delta)
  Wh <- sum(Z*sapply(1:n, function(x) phi(x,0.5+i,T)))
inc.ratio <- abs(Wh-W)/Delta
plot(Delta,inc.ratio,type="l",log="y",xlab=expression(Delta*t),
     ylab=expression(abs(W(0.5+Delta*t)-W(0.5))/Delta*t))
max(inc.ratio,na.rm=T)
