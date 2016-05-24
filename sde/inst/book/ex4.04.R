require(sde)
# ex4.04.R
set.seed(123)
theta <- c(6,2,1)
X <- sde.sim(X0 = rsCIR(1, theta), model="CIR", theta=theta,N=1000,delta=0.1)

f <-function(x) dsCIR(x, theta)
b <- function(x) theta[1]-theta[2]*x
sigma <- function(x) theta[3]*sqrt(x)
  
minX <- min(X)
maxX <- max(X)

par(mfrow=c(2,1))
curve(b,minX,maxX,main="drift coefficient")
lines(ksdrift(X),lty=3,col="red",lwd=2)

curve(sigma,minX, maxX,main="diffusion coefficient")
lines(ksdiff(X),lty=3,col="red",lwd=2)

 
 
 
 
 