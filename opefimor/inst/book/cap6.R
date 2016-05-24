###################################################
### chunk number 1: 
###################################################
#line 6 "cap6.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: callput-plot
###################################################
#line 35 "cap6.Rnw"
f.call <- function(x) sapply(x, function(x) max(c(x-K, 0)))
f.put <- function(x) sapply(x, function(x) max(c(K-x, 0)))
K <- 1
curve(f.call, 0, 2, main="Payoff functions", col="blue",lty=1,lwd=1,ylab=expression(f(x)))
curve(f.put, 0, 2, col="black", add=TRUE, lty=2,lwd=2)
legend(0.9,0.8, c("call","put"), lty=c(1,2), col=c("blue", "black"),lwd=c(2,2))


###################################################
### chunk number 3: 
###################################################
#line 45 "cap6.Rnw"
#line 35 "cap6.Rnw#from line#45#"
f.call <- function(x) sapply(x, function(x) max(c(x-K, 0)))
f.put <- function(x) sapply(x, function(x) max(c(K-x, 0)))
K <- 1
curve(f.call, 0, 2, main="Payoff functions", col="blue",lty=1,lwd=1,ylab=expression(f(x)))
curve(f.put, 0, 2, col="black", add=TRUE, lty=2,lwd=2)
legend(0.9,0.8, c("call","put"), lty=c(1,2), col=c("blue", "black"),lwd=c(2,2))
#line 46 "cap6.Rnw"


###################################################
### chunk number 4: callprice
###################################################
#line 730 "cap6.Rnw"
call.price <- function(x=1, t=0, T=1, r=1, sigma=1, K=1){	
	d2 <- (log(x/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
	d1 <- d2 + sigma*sqrt(T-t)
	x*pnorm(d1) - K * exp(-r*(T-t)) * pnorm(d2)
}


###################################################
### chunk number 5: putprice
###################################################
#line 738 "cap6.Rnw"
put.price <- function(x=1, t=0, T=1, r=1, sigma=1, K=1){	
	d2 <- (log(x/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
	d1 <- d2 + sigma*sqrt(T-t)
	K * exp(-r*(T-t)) * pnorm(-d2) - x*pnorm(-d1) 
}


###################################################
### chunk number 6: 
###################################################
#line 746 "cap6.Rnw"
S0 <- 100
K <- 110
r <- 0.05
T <- 1/4
sigma <- 0.25
C <- call.price(x=S0, t=0, T=T, r=r, K=K, sigma=sigma) 
C


###################################################
### chunk number 7: 
###################################################
#line 756 "cap6.Rnw"
P <- put.price(x=S0, t=0, T=T, r=r, K=K, sigma=sigma) 
P


###################################################
### chunk number 8: 
###################################################
#line 761 "cap6.Rnw"
C - S0+K*exp(-r*T)


###################################################
### chunk number 9: fOptions
###################################################
#line 765 "cap6.Rnw"
require(fOptions)


###################################################
### chunk number 10: 
###################################################
#line 769 "cap6.Rnw"
GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)


###################################################
### chunk number 11: 
###################################################
#line 773 "cap6.Rnw"
GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price


###################################################
### chunk number 12: 
###################################################
#line 778 "cap6.Rnw"
GBSOption(TypeFlag = "p", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price


###################################################
### chunk number 13: 
###################################################
#line 792 "cap6.Rnw"
a <- 5
a*GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
GBSOption(TypeFlag = "c", S = a*S0, X = a*K, Time = T, r = r, b = r, sigma = sigma)@price


###################################################
### chunk number 14: 
###################################################
#line 808 "cap6.Rnw"
MCPrice <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f){
	h <- function(m){  
            u <- rnorm(m/2)
            tmp <- c(x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
	        x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u)))
            mean( sapply(tmp, function(xx) f(xx)) )
        }
        p <- h(M)
        p*exp(-r*(T-t))		   
}


###################################################
### chunk number 15: 
###################################################
#line 821 "cap6.Rnw"
S0 <- 100
K <- 110
r <- 0.05
T <- 1/4
sigma <- 0.25
GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
f <- function(x) max(0, x-K)
set.seed(123)
M <- 1000
MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)
set.seed(123)
M <- 50000
MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)
set.seed(123)
M <- 1000000
MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)


###################################################
### chunk number 16: 
###################################################
#line 843 "cap6.Rnw"
MCPrice <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f){
       require(foreach)
       h <- function(m){  
            u <- rnorm(m/2)
            tmp <- c(x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
	        x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u)))
            mean( sapply(tmp, function(xx) f(xx)) )
        }
        nodes <- getDoParWorkers()
        p <- foreach(m = rep(M/nodes, nodes), .combine="c") %dopar% h(m)
        p <- mean(p)
        p*exp(-r*(T-t))		   
}


###################################################
### chunk number 17: 
###################################################
#line 860 "cap6.Rnw"
f <- function(x) max(0, x-110)
set.seed(123)
M <- 50000
MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)


###################################################
### chunk number 18: 
###################################################
#line 869 "cap6.Rnw"
require(snowfall)
sfInit(parallel = TRUE, cpus = 2)
cl <- sfGetCluster()
clusterSetupRNG(cl, seed = rep(123, 2))


###################################################
### chunk number 19: 
###################################################
#line 876 "cap6.Rnw"
require(foreach)
require(doSNOW)
registerDoSNOW(cl)


###################################################
### chunk number 20: 
###################################################
#line 882 "cap6.Rnw"
getDoParWorkers()


###################################################
### chunk number 21: 
###################################################
#line 886 "cap6.Rnw"
M <- 50000
MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)


###################################################
### chunk number 22: MCruns
###################################################
#line 896 "cap6.Rnw"
set.seed(123)
m <- c(10, 50, 100, 150, 200, 250, 500, 1000)
p1 <- NULL
err <- NULL
nM <- length(m)
repl <- 100
mat <- matrix(, repl, nM)
for(k in 1:nM){
 tmp <- numeric(repl)
 for(i in 1:repl)
  tmp[i] <- MCPrice(x=S0, t=0, T=T, r=r, sigma,  M=m[k], f=f)
 mat[,k] <- tmp 
 p1 <- c(p1, mean(tmp))
 err <- c(err, sd(tmp)) 
}
colnames(mat) <- m


###################################################
### chunk number 23: MCspeed eval=FALSE
###################################################
## #line 915 "cap6.Rnw"
## p0 <- GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
## minP <- min(p1-err)
## maxP <- max(p1+err)
## plot(m, p1, type="n", ylim=c(minP,maxP),axes=F,ylab="MC price",xlab="MC replications")
## lines(m, p1+ err, col="blue")
## lines(m, p1-err, col="blue")
## axis(2, p0, "B&S price")
## axis(1, m)
## boxplot(mat,add=TRUE,at=m, boxwex=15,col="orange",axes=F)
## points(m, p1, col="blue",lwd=3,lty=3)
## abline(h=p0, lty=2, col="red",lwd=3)


###################################################
### chunk number 24: MCspeedplot
###################################################
#line 930 "cap6.Rnw"
par(mar=c(4,3,0,0))
#line 915 "cap6.Rnw#from line#931#"
p0 <- GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
minP <- min(p1-err)
maxP <- max(p1+err)
plot(m, p1, type="n", ylim=c(minP,maxP),axes=F,ylab="MC price",xlab="MC replications")
lines(m, p1+ err, col="blue")
lines(m, p1-err, col="blue")
axis(2, p0, "B&S price")
axis(1, m)
boxplot(mat,add=TRUE,at=m, boxwex=15,col="orange",axes=F)
points(m, p1, col="blue",lwd=3,lty=3)
abline(h=p0, lty=2, col="red",lwd=3)
#line 932 "cap6.Rnw"


###################################################
### chunk number 25: 
###################################################
#line 940 "cap6.Rnw"
sfStop()
registerDoSEQ()


###################################################
### chunk number 26: sens2K eval=FALSE
###################################################
## #line 993 "cap6.Rnw"
## S0 <- 100
## r <- 0.01
## T <- 100
## p <- function(sigma) call.price(x=S0, t=0, T=T, r=r, K=K, sigma=sigma) 
## K <- 80
## curve(p, 0, 1,xlab=expression(sigma),ylab=expression(P[t])) 
## K <- 100
## curve(p, 0, 1, add=TRUE, lty=2) 
## K <- 150
## curve(p, 0, 1,add=TRUE, lty=3)
## legend(0.5, 90, c("K=80", "K=100", "K=150"), lty=1:3)


###################################################
### chunk number 27: sens2Kplot
###################################################
#line 1008 "cap6.Rnw"
par(mar=c(4,4,1,0))
#line 993 "cap6.Rnw#from line#1009#"
S0 <- 100
r <- 0.01
T <- 100
p <- function(sigma) call.price(x=S0, t=0, T=T, r=r, K=K, sigma=sigma) 
K <- 80
curve(p, 0, 1,xlab=expression(sigma),ylab=expression(P[t])) 
K <- 100
curve(p, 0, 1, add=TRUE, lty=2) 
K <- 150
curve(p, 0, 1,add=TRUE, lty=3)
legend(0.5, 90, c("K=80", "K=100", "K=150"), lty=1:3)
#line 1010 "cap6.Rnw"


###################################################
### chunk number 28: sens2T eval=FALSE
###################################################
## #line 1021 "cap6.Rnw"
## S0 <- 100
## r <- 0.01
## K <- 100
## T <- 10
## curve(p, 0, 1,xlab=expression(sigma),ylab=expression(P[t]),ylim=c(0,100)) 
## T <- 50
## curve(p, 0, 1, add=TRUE, lty=2) 
## T <- 100
## curve(p, 0, 1,add=TRUE, lty=3)
## legend(0.5, 40, c("T=10", "T=50", "T=100"), lty=1:3)


###################################################
### chunk number 29: sens2Tplot
###################################################
#line 1035 "cap6.Rnw"
par(mar=c(4,4,1,0))
#line 1021 "cap6.Rnw#from line#1036#"
S0 <- 100
r <- 0.01
K <- 100
T <- 10
curve(p, 0, 1,xlab=expression(sigma),ylab=expression(P[t]),ylim=c(0,100)) 
T <- 50
curve(p, 0, 1, add=TRUE, lty=2) 
T <- 100
curve(p, 0, 1,add=TRUE, lty=3)
legend(0.5, 40, c("T=10", "T=50", "T=100"), lty=1:3)
#line 1037 "cap6.Rnw"


###################################################
### chunk number 30: sensd2T eval=FALSE
###################################################
## #line 1155 "cap6.Rnw"
## delta <- function(x){
## 	d2 <- (log(x/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
## 	d1 <- d2 + sigma*sqrt(T-t)
## 	pnorm(d1)
## }
## 
## r <- 0.01
## K <- 100
## T <- 100
## sigma <- 0.05
## t <- 1
## curve( delta, 0, 200, xlab=expression(S[t]), ylab=expression(delta=a^H(t)))
## t <- 50
## curve( delta, 0, 200, lty=2,add=TRUE)
## t <- 99.5
## curve( delta, 0, 200,lty=3,add=TRUE)
## legend(150, 0.6, c("t=1", "t=50", "t=99.5"), lty=1:3)


###################################################
### chunk number 31: sensd2Tplot
###################################################
#line 1176 "cap6.Rnw"
par(mar=c(4,4,1,0))
#line 1155 "cap6.Rnw#from line#1177#"
delta <- function(x){
	d2 <- (log(x/K) + (r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t))
	d1 <- d2 + sigma*sqrt(T-t)
	pnorm(d1)
}

r <- 0.01
K <- 100
T <- 100
sigma <- 0.05
t <- 1
curve( delta, 0, 200, xlab=expression(S[t]), ylab=expression(delta=a^H(t)))
t <- 50
curve( delta, 0, 200, lty=2,add=TRUE)
t <- 99.5
curve( delta, 0, 200,lty=3,add=TRUE)
legend(150, 0.6, c("t=1", "t=50", "t=99.5"), lty=1:3)
#line 1178 "cap6.Rnw"


###################################################
### chunk number 32: 
###################################################
#line 1257 "cap6.Rnw"
r <- 0.01
K <- 100
T <- 100
sigma <- 0.05
t <- 10
St <- 70
h <- 1e-2
delta.num <- function(x) (call.price(x=x+h,t=t,T=T,sigma=sigma,r=r,K=K)
-call.price(x=x,t=t,T=T,sigma=sigma,r=r,K=K))/h
delta(St)
delta.num(St)


###################################################
### chunk number 33: 
###################################################
#line 1275 "cap6.Rnw"
delta.num2 <- function(x) (call.price(x=x+h,t=t,T=T,sigma=sigma,r=r,K=K)
 -call.price(x=x-h,t=t,T=T,sigma=sigma,r=r,K=K))/(2*h)

delta(St)
delta.num2(St)


###################################################
### chunk number 34: 
###################################################
#line 1283 "cap6.Rnw"
delta(St)
h <- 1e-3
delta.num(St)
h <- 1e-4
delta.num(St)
h <- 1e-5
delta.num(St)


###################################################
### chunk number 35: MCdelta
###################################################
#line 1305 "cap6.Rnw"
MCdelta <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f){
       h <-function(m){
    	 u <- rnorm(M/2)
	 tmp <-  c(x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
	  x*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u))) 
        g <- function(z) f(z) * (log(z/x) - (r-0.5*sigma^2)*(T-t))/(x*sigma^2*(T-t))
        mean( sapply(tmp, function(z) g(z)) )
        }
        nodes <- getDoParWorkers()
        p <- foreach(m = rep(M/nodes, nodes), .combine="c") %dopar% h(m)
        p <- mean(p)
        p*exp(-r*(T-t))		   
}


###################################################
### chunk number 36: MCdeltaTest
###################################################
#line 1321 "cap6.Rnw"
r <- 0.01
K <- 100
T <- 100
t <- 10
sigma <- 0.05
S0 <- 70
f <- function(x) max(0, x-100)
delta(S0)
set.seed(123)
M <- 10000
MCdelta(x=S0, t=0, T=T, r=r, sigma,  M=M, f=f)


###################################################
### chunk number 37: 
###################################################
#line 1337 "cap6.Rnw"
r <- 0.01
K <- 100
T <- 100
t <- 0
sigma <- 0.05
S0 <- 70
delta(S0)
h <- 1e-3
p1 <- MCPrice(x=S0+h, t=0, T=T, r=r, sigma,  M=M, f=f)
p2 <- MCPrice(x=S0-h, t=0, T=T, r=r, sigma,  M=M, f=f)
p1
p2
(p1-p2)/(2*h)


###################################################
### chunk number 38: 
###################################################
#line 1354 "cap6.Rnw"
delta(S0)
set.seed(123)
p1 <- MCPrice(x=S0+h, t=0, T=T, r=r, sigma,  M=M, f=f)
set.seed(123)
p2 <- MCPrice(x=S0-h, t=0, T=T, r=r, sigma,  M=M, f=f)
p1
p2
(p1-p2)/(2*h)


###################################################
### chunk number 39: MCdelta2
###################################################
#line 1365 "cap6.Rnw"
MCdelta2 <- function(x=1, t=0, T=1, r=1, sigma=1, M=1000, f,dx=1e-3){
       h <-function(m){
    	 u <- rnorm(M/2)
	 tmp1 <-  c((x+dx)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
	  (x+dx)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u))) 
	 tmp2 <-  c((x-dx)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*u),
	  (x-dx)*exp((r-0.5*sigma^2)*(T-t)+sigma*sqrt(T-t)*(-u))) 
          mean( sapply(tmp1, function(x) f(x)) - sapply(tmp2, function(x) f(x)) )/(2*dx)
        }
        nodes <- getDoParWorkers()
        p <- foreach(m = rep(M/nodes, nodes), .combine="c") %dopar% h(m)
        p <- mean(p)
        p*exp(-r*(T-t))		   
}


###################################################
### chunk number 40: 
###################################################
#line 1382 "cap6.Rnw"
set.seed(123)
MCdelta2(x=S0,t=0,T=T,r=r,sigma=sigma,f=f,M=M, dx=h)
delta(S0)


###################################################
### chunk number 41: 
###################################################
#line 1452 "cap6.Rnw"
r <- 0.01
K <- 100
T <- 100
t <- 10
sigma <- 0.05
S0 <- 70
GBSCharacteristics(TypeFlag = "c", S = S0, X = K, Time = T-t, r = r, b = r, sigma = sigma)


###################################################
### chunk number 42: 
###################################################
#line 1462 "cap6.Rnw"
GBSCharacteristics(TypeFlag = "p", S = S0, X = K, Time = T-t, r = r, b = r, sigma = sigma)


###################################################
### chunk number 43: 
###################################################
#line 1772 "cap6.Rnw"
MCAsian <- function(S0=100, K=100, t=0, T=1, mu=0.1, sigma=0.1, r=0.1, N=100, M=1000){
   require(foreach)
   h <- function(x){
      require(sde)
      z <- colMeans(sde.sim(X0=S0, model="BS", theta=c(mu,sigma), M=x, N=N))
      f <- function(x) max(x-K,0)
      p0 <- mean(sapply(z, f))
   }
   nodes <- getDoParWorkers()
   p <- foreach(m = rep(M/nodes, nodes), .combine="c") %dopar% h(m)
   p <- mean(p)
   p*exp(-r*(T-t))
}


###################################################
### chunk number 44: 
###################################################
#line 1788 "cap6.Rnw"
M <- 5000
mu <- 0.1
s <- 0.5
K <- 110
r <- 0.01
T <- 1
S0 <- 100
set.seed(123)
p0 <- MCAsian(S0=S0, K=K, t=0, T=T, mu=mu, sigma=s, r=r, N=250, M=M)
p0


###################################################
### chunk number 45: 
###################################################
#line 1889 "cap6.Rnw"
require(yuima)


###################################################
### chunk number 46: 
###################################################
#line 1892 "cap6.Rnw"
AEAsian <- function(S0=100, K=100, t=0, T=1, mu=0.1, e=0.1, r=0.1, N=1000){
   require(yuima)
   diff.matrix <- matrix( c("x*e"), 1,1)
   model <- setModel(drift = c("mu*x"), diffusion = diff.matrix)
   xinit <- S0
   f <- list( expression(x/T), expression(0))
   F <- 0
   yuima <- setYuima(model = model, sampling = setSampling(Terminal=T, n=1000))
   yuima <- setFunctional( yuima, f=f,F=F, xinit=xinit,e=e)
   F0 <- F0(yuima)
   rho <- expression(0)
   get_ge <- function(x,epsilon,K,F0){
    tmp <- (F0 - K) + (epsilon * x) 
    tmp[(epsilon * x) < (K-F0)] <- 0
    return( tmp )
   }
   epsilon <- e  # noise level
   g <- function(x) {
    tmp <- (F0 - K) + (epsilon * x) 
    tmp[(epsilon * x) < (K-F0)] <- 0
    tmp
   }
   
   asymp <- asymptotic_term(yuima, block=10, rho, g)
   exp(-r*T)*(asymp$d0 + e * asymp$d1 )
}


###################################################
### chunk number 47: 
###################################################
#line 1922 "cap6.Rnw"
p1 <- AEAsian(S0=S0, K=K, t=0, T=T, mu=mu, e=s, r=r)
p1


###################################################
### chunk number 48: 
###################################################
#line 1931 "cap6.Rnw"
require(fExoticOptions)
p2 <- TurnbullWakemanAsianApproxOption("c", S = S0, SA = S0, X = K,
       Time = T, time = T, tau = 0.0, r = r, b = r, sigma = s)@price
p2


###################################################
### chunk number 49: 
###################################################
#line 1938 "cap6.Rnw"
require(fAsianOptions)
p3 <- GemanYorAsianOption("c", S = S0, X = K, Time = T,  r = r, sigma = s, doprint = FALSE)$price
p3
p4 <- VecerAsianOption("c", S = S0, X = K, Time = T, r = r, sigma = s, table = NA, 
 nint = 800, eps = 1.0e-8, dt = 1.0e-10)
p4
p5 <- ZhangAsianOption("c", S = S0, X = K, Time = T, r = r, sigma = s, table = NA, correction = TRUE, nint = 800,  eps = 1.0e-8, dt = 1.0e-10)
p5


###################################################
### chunk number 50: 
###################################################
#line 1973 "cap6.Rnw"
require(fImport)
S <- yahooSeries("ATL.MI", from="2004-07-23", to="2005-05-13")
head(S)
Close <- S[,"ATL.MI.Close"]


###################################################
### chunk number 51: figATL eval=FALSE
###################################################
## #line 1980 "cap6.Rnw"
## require(quantmod)
## chartSeries(Close, theme="white")


###################################################
### chunk number 52: 
###################################################
#line 1986 "cap6.Rnw"
#line 1980 "cap6.Rnw#from line#1986#"
require(quantmod)
chartSeries(Close, theme="white")
#line 1987 "cap6.Rnw"


###################################################
### chunk number 53: 
###################################################
#line 1994 "cap6.Rnw"
X <- returns(Close)
Delta <- 1/252
sigma.hat <- sqrt(var(X)/Delta)[1,1]
sigma.hat


###################################################
### chunk number 54: 
###################################################
#line 2004 "cap6.Rnw"
S0 <- Close[1]
K <- 23
T <- 15*Delta
r <- 0.02074
sigma.hat <- as.numeric(sigma.hat)
require(fOptions)
p0 <- GBSOption("c", S=S0, X=K, Time=T, r=r, b=r, sigma=sigma.hat)@price
p0


###################################################
### chunk number 55: 
###################################################
#line 2022 "cap6.Rnw"
p <- 0.0004
sigma.imp <- GBSVolatility(p, "c", S=S0, X=K, Time=T, r=r, b=r)
sigma.imp


###################################################
### chunk number 56: figATL2 eval=FALSE
###################################################
## #line 2030 "cap6.Rnw"
## require(sde)
## cp <- cpoint(as.ts(X))
## cp
## time(X)[cp$k0]
## plot(X)
## abline(v=time(X)[cp$k0], lty=3)


###################################################
### chunk number 57: 
###################################################
#line 2043 "cap6.Rnw"
par(mar=c(3,2,0,0))
#line 2030 "cap6.Rnw#from line#2044#"
require(sde)
cp <- cpoint(as.ts(X))
cp
time(X)[cp$k0]
plot(X)
abline(v=time(X)[cp$k0], lty=3)
#line 2045 "cap6.Rnw"


###################################################
### chunk number 58: smile
###################################################
#line 2073 "cap6.Rnw"
S <- yahooSeries("AAPL", from="2009-01-02", to="2009-04-23")
Close <- S[,"AAPL.Close"]
X <- returns(Close)
Delta <- 1/252
sigma.hat <- sqrt(var(X)/Delta)
sigma.hat

Pt <- c(22.2, 18.4, 15.02, 11.9, 9.2, 7, 5.2, 3.6, 2.62, 1.76, 1.28, 
0.8, 0.53, 0.34, 0.23, 0.15, 0.09, 0.1)

K <- c(105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 
165, 170, 175, 180, 185, 190)

S0 <- 123.90
nP <- length(Pt)
T <- 60*Delta
r <- 0.056

smile <- sapply(1:nP, function(i) GBSVolatility(Pt[i], "c", S=S0, X=K[i], Time=T, r=r, b=r))


###################################################
### chunk number 59: figSMILE eval=FALSE
###################################################
## #line 2095 "cap6.Rnw"
## vals <- c(smile, sigma.hat)
## 
## plot(K,smile,type="l",ylim=c(min(vals,na.rm=TRUE), max(vals,na.rm=TRUE)),main="")
## abline(v=S0,lty=3,col="blue")
## abline(h=sigma.hat,lty=3,col="red")
## axis(2, sigma.hat, expression(hat(sigma)), col="red")


###################################################
### chunk number 60: 
###################################################
#line 2105 "cap6.Rnw"
par(mar=c(5,3,0,0))
#line 2095 "cap6.Rnw#from line#2106#"
vals <- c(smile, sigma.hat)

plot(K,smile,type="l",ylim=c(min(vals,na.rm=TRUE), max(vals,na.rm=TRUE)),main="")
abline(v=S0,lty=3,col="blue")
abline(h=sigma.hat,lty=3,col="red")
axis(2, sigma.hat, expression(hat(sigma)), col="red")
#line 2107 "cap6.Rnw"


