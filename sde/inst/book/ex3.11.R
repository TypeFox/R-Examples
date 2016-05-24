# ex3.11.R
require(sde)
require(stats4)
# ex3.11.R
K.est <- function(x) {
  n.obs <- length(x)
  n.obs/(2*(sum(x^2)))
}

LS.est <- function(x) {
  n <- length(x) -1
  k.sum <- sum(x[1:n]*x[2:(n+1)])
  dt <- deltat(x)
  ifelse(k.sum>0, -log(k.sum/sum(x[1:n]^2))/dt, NA)
}

MLE.est <- function(y, lower=0, upper=Inf){
 n <- length(y) - 1
 Dt <- deltat(y)
 Y <- y[2:(n+1)]
 g <- function(theta){
  ss <- sqrt((1-exp(-2*Dt*theta))/(2*theta))
  X <- y[1:n]*exp(-theta*Dt)
  lik <- dnorm(Y,mean=X,sd=ss)  
  -sum(log(lik))
 } 
 tmp <- try(optim(runif(1), g, method="L-BFGS-B", lower=lower, upper=upper)$par) 
 if(class(tmp)=="try-error") tmp <- NA
 tmp
}


# ex3.11.R (cont)
theta0 <- 1
d <- expression( -1*x ) 
s <-  expression( 1 )

K <- 1000 # 1000 MC replications

set.seed(123)
kessler <- matrix(NA,K,9)
mle <- matrix(NA,K,9)
simple <- matrix(NA,K,9)
 x0 <- rnorm(K,sd=sqrt(1/(2*theta0)))
 sde.sim(X0=x0, drift=d, sigma=s, N=50000, delta=0.1, M=K)->X
for(k in 1:K){
 cat(".") 
 m <- 0
 for(Delta in c(0.4,1,5)){
  m <- m+1
  j <- 0
  for(n in c(200,500,1000)){
   j <- j+1
   X.win <- window(X[,k], start=0, end=n*Delta, deltat=Delta)
   kessler[k,m+3*(j-1)] <- K.est(X.win)
   simple[k,m+3*(j-1)] <- LS.est(X.win)
   mle[k,m+3*(j-1)] <- MLE.est(X.win)
  }
 }
 cat(sprintf(" %3.3d / %3.3d completed\n",k,K))
}

# ex3.11.R (cont)
S1 <- apply(simple,2,function(x) mean(x,na.rm=T))
K1 <- apply(kessler,2,function(x) mean(x,na.rm=T))
M1 <- apply(mle,2,function(x) mean(x,na.rm=T))
A <- cbind(S1,K1,M1)
matrix(as.numeric(sprintf("%3.2f",A)),9,3)

S2 <- apply(simple,2,function(x) sd(x,na.rm=T))
K2 <- apply(kessler,2,function(x) sd(x,na.rm=T))
M2 <- apply(mle,2,function(x) sd(x,na.rm=T))
B <- cbind(S2,K2,M2)
matrix(as.numeric(sprintf("%3.2f",B)),9,3)

Delta <- c(0.4,1,5)
v0 <- 2*theta0^2 * (1+exp(-2*theta0*Delta))/(1-
  exp(-2*theta0*Delta))

sprintf("%3.2f",sqrt(v0[1]/c(200,500,1000)))
sprintf("%3.2f",sqrt(v0[2]/c(200,500,1000)))
sprintf("%3.2f",sqrt(v0[3]/c(200,500,1000)))

# valid cases for the LS estimator
apply(simple, 2, function(x) length(which(!is.na(x))))




require(sde)
# ex3.11.R (cont)
# CIR-Model
# theta = -1
# alpha = 10
# sigma = 1
# x0 = 10
set.seed(123); 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1000,delta=0.1) -> X

# estimator for alpha
(sum(X)^2)/(2*(length(X)*sum(X^2)-sum(X)^2))

# estimator for theta
(-length(X)*sum(X))/(2*(length(X)*sum(X^2)-sum(X)^2))


