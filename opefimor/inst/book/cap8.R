###################################################
### chunk number 1: 
###################################################
#line 6 "cap9.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: 
###################################################
#line 281 "cap9.Rnw"
FFTcall.price <- function(phi, S0, K, r, T, alpha=1, N=2^12, eta=0.25 ){
 m <- r-log(phi(-(1i)))  # mean-correcting mart. measure
 phi.tilde <- function(u) (phi(u)*exp((1i)*u*m))^T # EMM cf
 psi <- function(v) exp(-r*T)*phi.tilde((v-(alpha+1)*(1i)))/(alpha^2+alpha- v^2+(1i)*(2*alpha+1)*v)

 lambda <- (2*pi)/(N*eta) #(23) 
 b <- 1/2*N*lambda	#(20) 
 ku <- -b+lambda*(0:(N-1)) #(19) 
 v <- eta*(0:(N-1))	#above (17) 
 tmp <- exp((1i)*b*v)*psi(v)*eta*(3+(-1)^(1:N)-((1:N)-1==0))/3 #(24) 
 ft <- fft(tmp)
 res <- exp(-alpha*ku)*ft/pi #(24) 
 inter <- spline(ku,Re(res),xout=log(K/S0)) #Interpolate to obtain the value 
 return(inter$y*S0) 
}


###################################################
### chunk number 3: 
###################################################
#line 304 "cap9.Rnw"
phiBS <- function(u) exp(1i*u*(mu-0.5*sigma^2) - 0.5*sigma^2*u^2)


###################################################
### chunk number 4: 
###################################################
#line 308 "cap9.Rnw"
S0 <- 100
K <- 110
r <- 0.05
T <- 1/4
sigma <- 0.25
mu <- 1

require(fOptions)

GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price
FFTcall.price(phiBS, S0=S0, K=K, r=r, T=T)


###################################################
### chunk number 5: fftPlot eval=FALSE
###################################################
## #line 322 "cap9.Rnw"
## K.seq <- seq(100, 120, length=100)
## exactP <- NULL
## fftP <- NULL
## for(K in K.seq){
##  exactP  <- c(exactP , GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price )
##  fftP  <- c(fftP , FFTcall.price(phiBS, S0=S0, K=K, r=r, T=T) )
## }
## plot(K.seq, exactP-fftP, type="l",xlab="strike price K")


###################################################
### chunk number 6: 
###################################################
#line 334 "cap9.Rnw"
par(mar=c(4,4,1,1))
#line 322 "cap9.Rnw#from line#335#"
K.seq <- seq(100, 120, length=100)
exactP <- NULL
fftP <- NULL
for(K in K.seq){
 exactP  <- c(exactP , GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma)@price )
 fftP  <- c(fftP , FFTcall.price(phiBS, S0=S0, K=K, r=r, T=T) )
}
plot(K.seq, exactP-fftP, type="l",xlab="strike price K")
#line 336 "cap9.Rnw"


###################################################
### chunk number 7: 
###################################################
#line 359 "cap9.Rnw"
theta <- -0.1436
nu <- 0.3
r <- 0.1
sigma <- 0.12136
T <- 1
K <- 101
S <- 100
alpha <- 1.65
phiVG <- function(u) {
 omega <- (1/nu)*(log(1-theta*nu-sigma^2*nu/2))
 tmp <- 1-1i*theta*nu*u + 0.5*sigma^2*u^2*nu
 tmp <- tmp^(-1/nu)
 exp(1i*u*log(S0)+u*(r+omega)*1i)*tmp
}


###################################################
### chunk number 8: 
###################################################
#line 376 "cap9.Rnw"
FFTcall.price(phiVG, S0=S0, K=K, r=r, T=T)


###################################################
### chunk number 9: 
###################################################
#line 385 "cap9.Rnw"
n <- 50000
t <- rgamma(n, shape=T/nu, scale=nu)
N <- rnorm(n, 0,1)
X <- theta*t + N*sigma*sqrt(t)
omega <- (1/nu)*(log(1-theta*nu-sigma^2*nu/2))
S <- S0*exp( r*T+omega*T+X )
payoff <- sapply(S, function(x) max(x-K,0))
mean(payoff)*exp(-r*T)


###################################################
### chunk number 10: fftPlot2 eval=FALSE
###################################################
## #line 396 "cap9.Rnw"
## K.seq <- seq(100, 120, length=100)
## mcP <- NULL
## fftP <- NULL
## for(K in K.seq){
##  t <- rgamma(n, shape=T/nu, scale=nu)
##  N <- rnorm(n, 0,1)
##  X <- theta*t + N*sigma*sqrt(t)
##  S <- S0*exp( r*T+omega*T+X )
##  payoffvec <- sapply(S, function(x) max(x-K,0))
##  tmp <- mean(payoffvec)*exp(-r*T)
## 
##  mcP  <- c(mcP , tmp )
##  fftP  <- c(fftP , FFTcall.price(phiVG, S0=S0, K=K, r=r, T=T) )
## }
## plot(K.seq, mcP-fftP, type="l",xlab="strike price K")


###################################################
### chunk number 11: 
###################################################
#line 417 "cap9.Rnw"
par(mar=c(4,4,1,1))
#line 396 "cap9.Rnw#from line#418#"
K.seq <- seq(100, 120, length=100)
mcP <- NULL
fftP <- NULL
for(K in K.seq){
 t <- rgamma(n, shape=T/nu, scale=nu)
 N <- rnorm(n, 0,1)
 X <- theta*t + N*sigma*sqrt(t)
 S <- S0*exp( r*T+omega*T+X )
 payoffvec <- sapply(S, function(x) max(x-K,0))
 tmp <- mean(payoffvec)*exp(-r*T)

 mcP  <- c(mcP , tmp )
 fftP  <- c(fftP , FFTcall.price(phiVG, S0=S0, K=K, r=r, T=T) )
}
plot(K.seq, mcP-fftP, type="l",xlab="strike price K")
#line 419 "cap9.Rnw"


###################################################
### chunk number 12: 
###################################################
#line 692 "cap9.Rnw"
V0 <- function(S0, K, T, r, s0, s1, L0, L1, p0){	
	m <- function(t) log(S0) +(d1-d0-0.5*(s0^2-s1^2))*t +(r-d1-0.5*s1^2)*T
	v <- function(t) (s0^2-s1^2)*t+s1^2*T	
	
	f <- function(y){
		y/(y+K) * (dnorm(log(y+K), m(p0*T), sqrt(v(p0*T)))*(1-exp(-L0*T))+
				   dnorm(log(y+K), m(T), sqrt(v(T)))*exp(-L0*T) )
	}	
	integrate(f, 0, Inf, subdivisions=1000)$value * exp(-r*T)
}

V1 <- function(S0, K, T, r, s0, s1, L0, L1, p0){	
	m <- function(t) log(S0) +(d1-d0-0.5*(s0^2-s1^2))*t +(r-d1-0.5*s1^2)*T
	v <- function(t) (s0^2-s1^2)*t+s1^2*T
	
	f <- function(y){
		y/(y+K) *( dnorm(log(y+K), m(p0*T), sqrt(v(p0*T)))*(1-exp(-L1*T))+
				  dnorm(log(y+K), m(0), sqrt(v(0)))*exp(-L1*T) )
	}	
	integrate(f, 0, Inf, subdivisions=1000)$value * exp(-r*T)
}


###################################################
### chunk number 13: MSWPlot eval=FALSE
###################################################
## #line 716 "cap9.Rnw"
## require(fOptions)
## r <- 0.1
## s0 <- 0.2
## s1 <- 0.4
## L0 <- 1
## L1 <- 1
## K <- 110
## S0 <- 100
## d1 <- 0
## d0 <- 0
## p0 <- 0.5 
## 
## tt <- seq(0,1, length=50)
## 
## pV0 <- NULL
## pV1 <- NULL
## pBS0 <- NULL
## pBS1 <- NULL
## for(T in tt){
## 	pV0 <- c(pV0, V0(S0, K, T, r, s0, s1, L0, L1, p0))
## 	pV1 <- c(pV1, V1(S0, K, T, r, s0, s1, L0, L1, p0))
## 	pBS0 <- c(pBS0, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
## 							  b = r, sigma = s0)@price)
## 	pBS1 <- c(pBS1, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
## 							  b = r, sigma = s1)@price)
## }
## matplot(tt, cbind(pV0,pV1,pBS0,pBS1),type="o", lty=rep(1,4), pch=c(0,1,15,16),cex=0.5,col=rep(1,4),main="",xlab=expression(T))
## legend(0.1,12, c(expression(tilde(V)[0]), expression(tilde(V)[1]), expression(BS[0]), expression(BS[1])), 
## lty=rep(1,4), pch=c(0,1,15,16),cex=0.75,col=rep(1,4))


###################################################
### chunk number 14: 
###################################################
#line 749 "cap9.Rnw"
par(mar=c(4.5,2.5,0.5,0.5))
#line 716 "cap9.Rnw#from line#750#"
require(fOptions)
r <- 0.1
s0 <- 0.2
s1 <- 0.4
L0 <- 1
L1 <- 1
K <- 110
S0 <- 100
d1 <- 0
d0 <- 0
p0 <- 0.5 

tt <- seq(0,1, length=50)

pV0 <- NULL
pV1 <- NULL
pBS0 <- NULL
pBS1 <- NULL
for(T in tt){
	pV0 <- c(pV0, V0(S0, K, T, r, s0, s1, L0, L1, p0))
	pV1 <- c(pV1, V1(S0, K, T, r, s0, s1, L0, L1, p0))
	pBS0 <- c(pBS0, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s0)@price)
	pBS1 <- c(pBS1, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s1)@price)
}
matplot(tt, cbind(pV0,pV1,pBS0,pBS1),type="o", lty=rep(1,4), pch=c(0,1,15,16),cex=0.5,col=rep(1,4),main="",xlab=expression(T))
legend(0.1,12, c(expression(tilde(V)[0]), expression(tilde(V)[1]), expression(BS[0]), expression(BS[1])), 
lty=rep(1,4), pch=c(0,1,15,16),cex=0.75,col=rep(1,4))
#line 751 "cap9.Rnw"


###################################################
### chunk number 15: 
###################################################
#line 762 "cap9.Rnw"
par(mar=c(4.5,2.5,0.5,0.5))
require(fOptions)
r <- 0.1
s0 <- 0.2
s1 <- 0.4
L0 <- 10
L1 <- 0.1
K <- 110
S0 <- 100
d1 <- 0
d0 <- 0
p0 <- 0.5 

tt <- seq(0,1, length=50)

pV0 <- NULL
pV1 <- NULL
pBS0 <- NULL
pBS1 <- NULL
for(T in tt){
	pV0 <- c(pV0, V0(S0, K, T, r, s0, s1, L0, L1, p0))
	pV1 <- c(pV1, V1(S0, K, T, r, s0, s1, L0, L1, p0))
	pBS0 <- c(pBS0, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s0)@price)
	pBS1 <- c(pBS1, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s1)@price)
}
matplot(tt, cbind(pV0,pV1,pBS0,pBS1),type="o", lty=rep(1,4), pch=c(0,1,15,16),cex=0.5,col=rep(1,4),main="",xlab=expression(T))
legend(0.1,12, c(expression(tilde(V)[0]), expression(tilde(V)[1]), expression(BS[0]), expression(BS[1])), 
lty=rep(1,4), pch=c(0,1,15,16),cex=0.75,col=rep(1,4))


###################################################
### chunk number 16: 
###################################################
#line 802 "cap9.Rnw"
par(mar=c(4.5,2.5,0.5,0.5))
require(fOptions)
r <- 0.1
s0 <- 0.2
s1 <- 0.4
L0 <- 10
L1 <- 10
K <- 110
S0 <- 100
d1 <- 0
d0 <- 0
p0 <- 0.5 

tt <- seq(0,1, length=50)

pV0 <- NULL
pV1 <- NULL
pBS0 <- NULL
pBS1 <- NULL
for(T in tt){
	pV0 <- c(pV0, V0(S0, K, T, r, s0, s1, L0, L1, p0))
	pV1 <- c(pV1, V1(S0, K, T, r, s0, s1, L0, L1, p0))
	pBS0 <- c(pBS0, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s0)@price)
	pBS1 <- c(pBS1, GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, 
							  b = r, sigma = s1)@price)
}
matplot(tt, cbind(pV0,pV1,pBS0,pBS1),type="o", lty=rep(1,4), pch=c(0,1,15,16),cex=0.5,col=rep(1,4),main="",xlab=expression(T))
legend(0.1,12, c(expression(tilde(V)[0]), expression(tilde(V)[1]), expression(BS[0]), expression(BS[1])), 
lty=rep(1,4), pch=c(0,1,15,16),cex=0.75,col=rep(1,4))


###################################################
### chunk number 17: 
###################################################
#line 866 "cap9.Rnw"
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
	
    ts(x, deltat=delta,start=0)	
}


###################################################
### chunk number 18: 
###################################################
#line 902 "cap9.Rnw"
T <- 0.75
r <- 0.05
s0 <- 0.2
s1 <- 0.4
L0 <- 1
L1 <- 1
K <- 80
S0 <- 85
d1 <- 0
d0 <- 0
p0 <- 0.5
Q <- matrix( c(-L0, L1, L0, -L1), 2, 2)


###################################################
### chunk number 19: 
###################################################
#line 917 "cap9.Rnw"
c(p0, p0) %*% Q


###################################################
### chunk number 20: 
###################################################
#line 921 "cap9.Rnw"
f <-function(x,a)  ifelse(a==0, (r-d0)*x, (r-d1)*x)
g <-function(x,a) ifelse(a==0, s0*x, s1*x)


###################################################
### chunk number 21: simMS0
###################################################
#line 931 "cap9.Rnw"
n <- 1000
set.seed(123)
nsim <- 1000
ST0 <- numeric(nsim)
for(i in 1:nsim){
 X <- simMSdiff(x0=S0, a0=0, S=0:1, delta=1/n, n=T*n, f, g, Q)
 ST0[i] <- X[length(X)]
}
v0 <- exp(-r*T)*pmax(ST0-K,0)
mean(v0)


###################################################
### chunk number 22: 
###################################################
#line 944 "cap9.Rnw"
mean(v0)
V0(S0, K, T, r, s0, s1, L0, L1, p0)


###################################################
### chunk number 23: simMS1
###################################################
#line 949 "cap9.Rnw"
set.seed(123)
ST1 <- numeric(nsim)
for(i in 1:nsim){
	X <- simMSdiff(x0=S0, a0=1, S=0:1, delta=1/n, n=T*n, f, g, Q)
	ST1[i] <- X[length(X)]
}
v1 <- exp(-r*T)*pmax(ST1-K,0)
mean(v1)
V1(S0, K, T, r, s0, s1, L0, L1, p0)


###################################################
### chunk number 24: 
###################################################
#line 994 "cap9.Rnw"
Vsmc <- function(S0, K, T, r, a0, S, Q, mu, sigma, M=1000){
    require(msm)
    C <- numeric(M)
    for(i in 1:M){
         MC <- sim.msm(Q, T, start=which(S==a0))
	 alpha <- S[MC$states]
	
	 t <- diff(MC$times)
	 s <- sapply(alpha[-length(alpha)], sigma)^2
	 m <- sapply(alpha[-length(alpha)], mu)
	 VT <- sum(t*s)
	 LT <- sum(t*m)
	
	 d1 <- (log(S0/K) + LT + 0.5*VT)/sqrt(VT)
	 d2 <- d1 - sqrt(VT)
	 RT <- r*T 
	 C[i] <- S0*exp(-(RT -LT))*pnorm(d1) - K*exp(-RT)*pnorm(d2)
	}
   mean(C)
}


###################################################
### chunk number 25: 
###################################################
#line 1017 "cap9.Rnw"
T <- 1
r <- 0.1
s0 <- 0.2
s1 <- 0.3
L0 <- 1
L1 <- 1
K <- 90
S0 <- 100
d1 <- 0
d0 <- 0
p0 <- 0.5
Q <- matrix( c(-L0, L1, L0, -L1), 2, 2)
sigma <- function(a) ifelse(a==0, s0, s1)
mu <- function(a) ifelse(a==0, r-d0, r-d1)


###################################################
### chunk number 26: 
###################################################
#line 1034 "cap9.Rnw"
Vsmc(S0, K, T, r, a0=0, 0:1, Q, mu, sigma)
V0(S0, K, T, r, s0, s1, L0, L1, p0)


###################################################
### chunk number 27: 
###################################################
#line 1039 "cap9.Rnw"
Vsmc(S0, K, T, r, a0=1, 0:1, Q, mu, sigma)
V1(S0, K, T, r, s0, s1, L0, L1, p0)


