require(sde)
# ex4.03.R
set.seed(123)
theta <- c(6,2,1)
X <- sde.sim(X0 = rsCIR(1, theta), model="CIR", theta=theta,N=1000,delta=1)
f <-function(x) dsCIR(x, theta)

h <- length(X)^(-1/5)*sd(X) 
K <-function(x) exp(-0.5*x^2)/sqrt(2*pi)
p <-function(x) sapply(x, function(x) mean(K((x-X)/h)))/h

curve(f,0,8,ylim=c(0,0.7))
curve(p,0,8, col="red", add=TRUE,lty=2)
lines(density(X,bw=h),col="green",lty=3)
  
# ex4.03.R (cont)
set.seed(123)
X <- sde.sim(X0 = rsCIR(1, theta), model="CIR", theta=theta,N=1000,delta=0.01)
h <- length(X)^(-1/5)*sd(X) 
curve(f,0,8,ylim=c(0,0.7))
curve(p,0,8, col="red", add=TRUE,lty=2)

# ex4.03.R (cont)
set.seed(123)
X <- sde.sim(X0 = rsCIR(1, theta), model="CIR", theta=theta,N=15000,delta=0.01)
h <- length(X)^(-1/5)*sd(X) 
curve(f,0,8,ylim=c(0,0.7))
curve(p,0,8, col="red", add=TRUE,lty=2)



