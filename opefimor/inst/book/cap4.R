###################################################
### chunk number 1: 
###################################################
#line 6 "cap4.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: 
###################################################
#line 76 "cap4.Rnw"
set.seed(123)
g <- function(x) x^2
a <- 0
b <- 2
c <- 0
d <- 4
A <- (b-a)*(d-c)
n <- 100000
x <- runif(n, a, b)
y <- runif(n, c, d)
A*sum(y < g(x))/n
integrate(g,a,b)


###################################################
### chunk number 3: 
###################################################
#line 93 "cap4.Rnw"
set.seed(123)
g <- function(x) x^2
a <- 0
b <- 2
c <- 0
d <- 4
A <- (b-a)*(d-c)
n <- 100000
x <- runif(n, a, b)
y <- runif(n, c, d)

par(mfrow=c(1,3))
par(mar=c(3,4,2,1))
for(m in c(1000, 10000, 100000)){
 which(y[1:m]<g(x[1:m])) -> idx
 curve(g,a,b,main=sprintf("n=%d, val=%.4f",m,A*length(idx)/m))
 points((x[1:m])[idx], (y[1:m])[idx], col="red",pch=".")
}


###################################################
### chunk number 4: 
###################################################
#line 146 "cap4.Rnw"
h <- 0.01
x0 <- 1
err <- 1 - ( (x0+h)^(x0+h) - x0^x0 )/h
err


###################################################
### chunk number 5: 
###################################################
#line 153 "cap4.Rnw"
h <- 0.001
x0 <- 1
err <- 1 - ( (x0+h)^(x0+h) - x0^x0 )/h
err


###################################################
### chunk number 6: 
###################################################
#line 178 "cap4.Rnw"
h <- 0.01
x0 <- 1
err <- 1 - ( (x0+h)^(x0+h) - (x0-h)^(x0-h) )/(2*h)
err
h^2
h <- 0.001
x0 <- 1
err <- 1 - ( (x0+h)^(x0+h) - (x0-h)^(x0-h) )/(2*h)
err
h^2


###################################################
### chunk number 7: 
###################################################
#line 221 "cap4.Rnw"
h <- 0.001
x0 <- 1
err <- 1 - ( 2 * (x0+h/2)^(x0+h/2) - (x0+h)^(x0+h) )
err


###################################################
### chunk number 8: 
###################################################
#line 228 "cap4.Rnw"
require(numDeriv)
f <- function(x) x^x
grad(f,x=1)
grad(f,x=1, method="simple")


###################################################
### chunk number 9: 
###################################################
#line 256 "cap4.Rnw"
polyroot(c(1, 2, 1))


###################################################
### chunk number 10:  eval=FALSE
###################################################
## #line 264 "cap4.Rnw"
## f <- function(x) 1+2*x+x^2
## uniroot(f,c(-2,2))


###################################################
### chunk number 11: 
###################################################
#line 269 "cap4.Rnw"
f <- function(x) 2*x+x^2 -2
uniroot(f,c(-2,2))


###################################################
### chunk number 12: 
###################################################
#line 283 "cap4.Rnw"
f <- function(x) 2*x+x^2 -2
nlm(f, 0)


###################################################
### chunk number 13: 
###################################################
#line 289 "cap4.Rnw"
f <- function (x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80


###################################################
### chunk number 14: 
###################################################
#line 294 "cap4.Rnw"
par(mar=c(3,3,1,1))
curve(f, -50,50,n=1000,main="")


###################################################
### chunk number 15: 
###################################################
#line 303 "cap4.Rnw"
nlm(f, 0)$estimate
nlm(f, 20)$estimate


###################################################
### chunk number 16: 
###################################################
#line 309 "cap4.Rnw"
res <- optim(50, f, method="SANN", control=list(maxit=20000, temp=20, parscale=20))
res
res <- optim(0, f, method="SANN", control=list(maxit=20000, temp=20, parscale=20))
res$par


###################################################
### chunk number 17: 
###################################################
#line 326 "cap4.Rnw"
set.seed(123)
lambda <- 0.8
T <- 10
avg <- lambda*T
avg
t <- 0
N <- 0
k <- 0
continue <- TRUE
while(continue){
 event <- rexp(1, lambda)
 if(sum(t) + event < T){
  k <- k +1
  N <- c(N,k)
  t <- c(t, event)
 } else {
  continue <- FALSE
  t <- cumsum(t)
  N <- c(N,k)
  t <- c(t,T)
 }
}
N
t


###################################################
### chunk number 18: ppfig
###################################################
#line 353 "cap4.Rnw"
set.seed(123)
t <- cumsum(c(0, rexp(10*avg, lambda)))
last <- which(t>T)[1]
t <- t[1:last]
t[last] <- T
N <- c(0, 1:(length(t)-2), length(t)-2)
N
t
plot(t, N,type="s", main="Poisson process",ylab=expression(N(t)),xlim=c(0,T))


###################################################
### chunk number 19: 
###################################################
#line 366 "cap4.Rnw"
#line 353 "cap4.Rnw#from line#366#"
set.seed(123)
t <- cumsum(c(0, rexp(10*avg, lambda)))
last <- which(t>T)[1]
t <- t[1:last]
t[last] <- T
N <- c(0, 1:(length(t)-2), length(t)-2)
N
t
plot(t, N,type="s", main="Poisson process",ylab=expression(N(t)),xlim=c(0,T))
#line 367 "cap4.Rnw"


###################################################
### chunk number 20: inppfig
###################################################
#line 425 "cap4.Rnw"
set.seed(123)
lambda <- 1.1
T <- 20
E <- 0
t <- 0
while(t<T){
  t <- t - 1/lambda * log(runif(1))
  if( runif(1) < sin(t)/lambda )
   E <- c(E, t) 
}
plot(E,0:(length(E)-1),type="s",ylim=c(-4,length(E)),ylab=expression(N(t)),xlab="t")
curve(-3+sin(x),0,20,add=TRUE,lty=2,lwd=2)


###################################################
### chunk number 21: 
###################################################
#line 441 "cap4.Rnw"
#line 425 "cap4.Rnw#from line#441#"
set.seed(123)
lambda <- 1.1
T <- 20
E <- 0
t <- 0
while(t<T){
  t <- t - 1/lambda * log(runif(1))
  if( runif(1) < sin(t)/lambda )
   E <- c(E, t) 
}
plot(E,0:(length(E)-1),type="s",ylim=c(-4,length(E)),ylab=expression(N(t)),xlab="t")
curve(-3+sin(x),0,20,add=TRUE,lty=2,lwd=2)
#line 442 "cap4.Rnw"


###################################################
### chunk number 22: 
###################################################
#line 458 "cap4.Rnw"
T <- 20
lambda <- 5
avg <- lambda*T
t <- cumsum(c(0, rexp(10*avg, lambda)))
last <- which(t>T)[1]
t <- t[1:last]
t[last] <- T
N <- c(0, 1:(length(t)-2),length(t)-2)
c <- 2
V0 <- sample(c(-c,+c),1)
ds <- diff(t)
nds <- length(ds)
x0 <- 0
X <- c(x0, x0+cumsum(V0*ds*(-1)^(1:nds)))


###################################################
### chunk number 23: telprocfig
###################################################
#line 475 "cap4.Rnw"
par(mfrow=c(2,1))
par(mar=c(3,4,0.5,0.1))
plot(t,X,type="l")
plot(t,N,type="s")


###################################################
### chunk number 24: 
###################################################
#line 484 "cap4.Rnw"
#line 475 "cap4.Rnw#from line#484#"
par(mfrow=c(2,1))
par(mar=c(3,4,0.5,0.1))
plot(t,X,type="l")
plot(t,N,type="s")
#line 485 "cap4.Rnw"


###################################################
### chunk number 25: 
###################################################
#line 492 "cap4.Rnw"
T <- 1
c <- 100
lambda <- c^2
avg <- lambda*T
t <- cumsum(c(0, rexp(10*avg, lambda)))
last <- which(t>T)[1]
t <- t[1:last]
t[last] <- T
N <- c(0, 1:(length(t)-2),length(t)-2)
V0 <- sample(c(-c,+c),1)
ds <- diff(t)
nds <- length(ds)
x0 <- 0
X <- c(x0, x0+cumsum(V0*ds*(-1)^(1:nds)))


###################################################
### chunk number 26: 
###################################################
#line 511 "cap4.Rnw"
par(mar=c(3,4,0.5,0.1))
plot(t,X,type="l")


###################################################
### chunk number 27: sdesim
###################################################
#line 550 "cap4.Rnw"
require(sde)
set.seed(123)
b <- expression(1-2*x)
s <- expression(sqrt(1+x^2)) 
X <- sde.sim(X0=2, T=1, drift=b, sigma=s)


###################################################
### chunk number 28: sdesim2
###################################################
#line 560 "cap4.Rnw"
X <- sde.sim(X0=2, T=1, drift=b, sigma=s, method="milstein")


###################################################
### chunk number 29: sdesim3
###################################################
#line 568 "cap4.Rnw"
X <- sde.sim(X0=2, T=1, model="BS", theta=c(1,0.5))


###################################################
### chunk number 30: sdesim4
###################################################
#line 572 "cap4.Rnw"
set.seed(123)
X <- sde.sim(X0=2, T=1, model="BS", theta=c(1,0.5), N=5000)


###################################################
### chunk number 31: sdesim4
###################################################
#line 576 "cap4.Rnw"
plot(X)


###################################################
### chunk number 32: 
###################################################
#line 581 "cap4.Rnw"
#line 576 "cap4.Rnw#from line#581#"
plot(X)
#line 582 "cap4.Rnw"


###################################################
### chunk number 33: setModel
###################################################
#line 659 "cap4.Rnw"
require(yuima)
sol <- c("x1","x2") # variable for numerical solution
b <- c("-3*x1","-x1-2*x2")   # drift vector 
s <- matrix(c("1","x1","0","3","x2","0"),2,3)  #  diffusion matrix
model <- setModel(drift = b, diffusion = s, solve.variable = sol)


###################################################
### chunk number 34: simulate
###################################################
#line 668 "cap4.Rnw"
set.seed(123)
X <- simulate(model, n=1000)


###################################################
### chunk number 35: yuimaplot1
###################################################
#line 673 "cap4.Rnw"
plot(X,plot.type="single", lty=c(1,3),ylab="X")


###################################################
### chunk number 36: 
###################################################
#line 678 "cap4.Rnw"
#line 673 "cap4.Rnw#from line#678#"
plot(X,plot.type="single", lty=c(1,3),ylab="X")
#line 679 "cap4.Rnw"


###################################################
### chunk number 37: 
###################################################
#line 769 "cap4.Rnw"
require(yuima)
modCP <- setModel(drift=c("-theta*x"), diffusion="sigma",
 jump.coeff="1", measure=list(intensity="10", df=list("dnorm(z, 0, 1)")),
 measure.type="CP", solve.variable="x")


###################################################
### chunk number 38: 
###################################################
#line 776 "cap4.Rnw"
set.seed(123)
X <- simulate(modCP, true.p=list(theta=1,sigma=3),n=1000)


###################################################
### chunk number 39: modCP eval=FALSE
###################################################
## #line 781 "cap4.Rnw"
## plot(X, main="I'm jumping!")


###################################################
### chunk number 40: 
###################################################
#line 787 "cap4.Rnw"
par(mar=c(2,2,2,0))
#line 781 "cap4.Rnw#from line#788#"
plot(X, main="I'm jumping!")
#line 789 "cap4.Rnw"


###################################################
### chunk number 41: 
###################################################
#line 798 "cap4.Rnw"
modPJ <- setModel(drift="-x", xinit=1, jump.coeff="1", 
  measure.type="code", measure=list(df="rIG(z, 1, 0.1)"))
set.seed(123)
Y <-  simulate(modPJ, Terminal=10, n=10000)


###################################################
### chunk number 42: 
###################################################
#line 807 "cap4.Rnw"
par(mar=c(2,2,2,0))
plot(Y, main="I'm jumping as well!")


###################################################
### chunk number 43: 
###################################################
#line 847 "cap4.Rnw"
simMarkov <- function(x0, n, x, P){
	mk <- numeric(n+1)
	mk[1] <- x0
	state <- which(x==x0)
	for(i in 1:n){
		mk[i+1] <- sample(x,1,prob=P[state,])
		state <- which(x==mk[i+1])
	}
	return(ts(mk))
}


###################################################
### chunk number 44: 
###################################################
#line 869 "cap4.Rnw"
P <- matrix(c(0.1, 0.5, 0.3, 0.1, 0.2, 0.3, 0.8, 0.3, 0.4), 3, 3) 
set.seed(123)
X <- simMarkov(1, 10, 1:3, P)
plot(X, type="s")


###################################################
### chunk number 45: 
###################################################
#line 878 "cap4.Rnw"
par(mar=c(4,4,1,1))
plot(X, type="s")


###################################################
### chunk number 46: 
###################################################
#line 888 "cap4.Rnw"
simMSdiff <- function(x0, a0, S, delta, n, f, g, Q){
   require(msm)
    P <- MatrixExp(delta*Q)
    alpha <- simMarkov(a0, n, S, P)
	x <- numeric(n+1)
	x[1] <- x0
	for(i in 1:n){
		A <- f(x[i], alpha[i])*delta
		B <- g(x[i], alpha[i])*sqrt(delta)*rnorm(1)
		x[i+1] <- x[i] + A + B 		
	}
	
    ts(x, deltat=delta, start=0)	
}


###################################################
### chunk number 47: 
###################################################
#line 928 "cap4.Rnw"
Q <- matrix( c(-6, 12, 6, -12), 2, 2)
n <- 1000

s0 <- 0.5
s1 <- 0.2

f <-function(x,a)  ifelse(a==0, (1+sin(x))*x, (1-cos(x))*x)
g <-function(x,a) ifelse(a==0, s0*x, s1*x)

set.seed(123)
X <- simMSdiff (x0=5, a0=0, S=0:1, delta=1/n, n=1000, f, g, Q)
plot(X)


###################################################
### chunk number 48: 
###################################################
#line 945 "cap4.Rnw"
par(mar=c(4,4,1,1))
plot(X)


###################################################
### chunk number 49: 
###################################################
#line 958 "cap4.Rnw"
set.seed(123)
g <- function(x) sqrt(1-x^2)
a <- -1
b <- +1
c <- 0
d <- 2
A <- (b-a)*(d-c)
n <- 100000
x <- runif(n, a, b)
y <- runif(n, c, d)
2*A*sum(y < g(x))/n
2*integrate(g,a,b)$value
pi


