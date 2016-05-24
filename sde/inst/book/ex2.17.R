# ex2.17.R
require(sde)

drift <- expression((3-x))
sigma <- expression(1.2*sqrt(x))
a <- 1.7
b <- 0.5
set.seed(123)
Y1 <- sde.sim(X0=a, drift=drift, sigma=sigma, T=1, delta=0.01)
Y2 <- sde.sim(X0=b, drift=drift, sigma=sigma, T=1, delta=0.01)
Y3 <- ts(rev(Y2), start=start(Y2), end=end(Y2),deltat=deltat(Y2))


id1 <- Inf
if(Y1[1]>=Y3[1]){
 if(!all(Y1>Y3))
  min(which(Y1 <= Y3))-1 -> id1
} else {
 if(!all(Y1<Y3))
 min(which(Y1 >= Y3))-1 -> id1
}
if(id1==0 || id1==length(Y1)) id1 <- Inf

par(mar=c(3,3,1,1))
par(mfrow=c(2,1)) 
plot(Y1, ylim=c(min(Y1,Y2), max(Y1,Y2)),col="green",lty=2)
lines(Y3,col="blue",lty=3)

if(id1==Inf ){
 cat("no crossing")
} else {
 plot(Y1, ylim=c(min(Y1,Y2), max(Y1,Y2)),col="green",lty=2)
lines(Y3,col="blue",lty=3)
B <- ts(c(Y1[1:id1], Y3[-(1:id1)]), start=start(Y1),end=end(Y1),frequency=frequency(Y1))
lines(B,col="red",lwd=2)
}

# ex2.17.R (cont.)
d <- expression((3-x))
s <- expression(1.2*sqrt(x))
par(mar=c(3,3,1,1))
par(mfrow=c(2,1)) 
set.seed(123)
X <- DBridge(x=1.7,y=0.5, delta=0.01, drift=d, sigma=s)
plot(X)
X <- DBridge(x=1,y=5, delta=0.01, drift=d, sigma=s)
plot(X)
